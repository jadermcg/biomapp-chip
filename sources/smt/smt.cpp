#include "smt.h"
#include "smt_utils.h"

std::condition_variable cv, wrt;
std::mutex mtx;
std::queue<const arma::SpMat<uint64_t>*> save_queue;
std::ofstream smtdb;
std::atomic<bool> done_processing(false);
std::atomic<bool> done_writing(false);

//'Creates batches.
//'@name computeBatchSizes.
//'@param n Number of sequences.
//'@param k Size ok kmers.
//'@param t Width of sequences.
//'@return A std::tuple of bsize (batch size), nb (number of batches and r (size of remainder)).
std::tuple<int, int> computeBatchSize(int n, int bsize) {
  int nb {n / bsize};
  int r  {n % bsize};
  
  return std::make_tuple(nb, r);
}

//'Use TBB Concurrent HashMap to create a hash for fasta.
//'@name createTBBhash_.
//'@param fasta Dataset of sequences.
//'@param hash Tbb Concurrent HashMap.
//'@param k Size of the kmer.
//'@param start Initial sequence in actual batch.
//'@param end Final sequence in actual batch.
void createTBBhash_(const std::vector<std::string> &fasta, tbb::concurrent_hash_map<std::string, uint64_t> &hash, const int k, const int start, const int end) {
  
 const int n = end - start;
 const int t = fasta[0].size();
 const int m = t - k + 1;
 
 uint64_t nr = m * n * k;
 
 for (size_t i = start; i < end; ++i) {
  const auto &seq = fasta[i];
  for (size_t j = 0; j < m; ++j) {
    const auto &kmer = seq.substr(j, k);
    tbb::concurrent_hash_map<std::string, uint64_t>::accessor acc;
    hash.insert(acc, kmer);
    acc->second += 1;
  } 
 }
}

//'Use TBB Concurrent HashMap to create a hash for fasta.
//'@name createTBBhash.
//'@param fasta Dataset of sequences.
//'@param k Size of the kmer.
//' @param bsize Size of the bacthes.
//'@return Tbb Concurrent HashMap with kmers and yours counts with size k.
tbb::concurrent_hash_map<std::string, uint64_t> createTBBhash(std::vector<std::string> &fasta, const int k, const int bsize) {
  int n = fasta.size();
  int t = fasta[0].size();

  tbb::concurrent_hash_map<std::string, uint64_t> hash;
  
  // Number of batches
  int nb, r;
  std::tie(nb, r) = computeBatchSize(n, bsize);
  
  //setupBuffer
  int ret = std::system("rm -Rf smt_data");
  ret = std::system("mkdir smt_data"); 

  tbb::parallel_for(tbb::blocked_range<size_t>(0, nb), [&](tbb::blocked_range<size_t> r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      int start = i * bsize;
      int end = (i + 1) * bsize;
      createTBBhash_(fasta, hash, k, start, end);
    }
  });

  return hash;
}

//'Creates Dense matrix from the sequences by Armadillo.
//'@name createDenseMT
//'@param fasta The Dataset of sequences.
//'@param k The Size of kmers.
//'@param start The initial sequence will be processed.
//'@param end The final sequences will be processed.
//'@return A arma::Mat<uint> with all kmers stored.
arma::Mat<uint64_t>* createDenseMT(const std::vector<std::string> &fasta, const int k, const int start, const int end) {
 const auto n {end - start};
 const auto t { fasta[0].size() };
 const auto m { t - k + 1 };
 
 auto nr { m * n * k };
 auto *MT = new arma::Mat<uint64_t>(nr, 6);
 uint64_t next {0};
 uint64_t current {0};
 uint64_t symbol {0};
 
 for (auto i {start} ; i < end; ++i) {
   auto &seq {fasta[i]};
   
   for (int j = 0; j < m; ++j) {
     uint64_t node {0};
     uint64_t index {0};
     
     for (auto l {0}; l < k; ++l) {
       auto c {seq[j + l]};

        switch(c) {
          case 'A': symbol = 0; break;
          case 'C': symbol = 1; break;
          case 'G': symbol = 2; break;
          case 'T': symbol = 3; break;
        }

       index = index * 4 + symbol;
       current = MT->at(node, symbol);
       
       if (current == 0) {
         next += 1;
         MT->at(node,symbol) = next;
         node = next;
       }
       
       else {
         node = current;
       }
     }
     
     MT->at(node, 4) += 1;
     MT->at(node, 5) = index;
   }
 }
 
 return MT;
}

//'Creates Dense matrix from the sequences by Armadillo with direct memory allocation.
//'@name createDenseMT2
//'@param fasta The Dataset of sequences.
//'@param k The Size of kmers.
//'@param start The initial sequence will be processed.
//'@param end The final sequences will be processed.
//'@return A arma::Mat<uint> with all kmers stored.
uint64_t* createDenseMT2(const std::vector<std::string> &fasta, const int k, const int start, const int end) {
 const auto n {end - start};
 const auto t {fasta[0].size()};
 const auto m  {t - k + 1};
 
 auto nr { m * n * k};
 auto *MT {new uint64_t [nr * 6]};
 auto next {0};
 auto current {0};
 auto symbol {0};
 
 for (auto i {start}; i < end; ++i) {
   const auto &seq {fasta[i]};
   
   for (auto j {0}; j < m; ++j) {
     uint64_t node = 0;
     uint64_t index = 0;
     
     for (auto l {0}; l < k; ++l) {
       char c = seq[j + l];

        switch(c) {
          case 'A': symbol = 0; break;
          case 'C': symbol = 1; break;
          case 'G': symbol = 2; break;
          case 'T': symbol = 3; break;
        }

       index = index * 4 + symbol;
       current = MT[node*6 + symbol];
       
       if (current == 0) {
         next += 1;
         MT[node*6 + symbol] = next;
         node = next;
       }
       
       else {
         node = current;
       }
     }
     
     MT[node*6 + 4] += 1;
     MT[node*6 + 5] = index;
   }
 }
 
 return MT;
}

//'Thread that saves the data in SMT.db.
//'@name saveMT.
void saveMT() {

   while (!done_processing.load() || !save_queue.empty()) {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [&]{ return !save_queue.empty() || done_processing.load(); });
    if (!save_queue.empty()) {
      auto *S = save_queue.front();
      save_queue.pop();
      lock.unlock();
      S->save(smtdb, arma::arma_binary);
      delete S;
    }
  }

  done_writing.store(true);
  wrt.notify_one();
}

//'Creates SMT matrix from the sequences.
//'@name createSparseMT.
//'@param fasta The Dataset of sequences.
//'@param k The Size of kmers.
void processMT(const std::vector<std::string> &fasta, const int k, const int bsize) {
 
  const auto n { fasta.size() };
  const auto t { fasta[0].size() };
  const auto m { t - k + 1 };
  
  //setupBuffer
  int ret { std::system("rm -Rf smt_data") };
  ret = std::system("mkdir smt_data");
  smtdb.open("smt_data/SMT.db", std::ios::binary);

  // Number of batches
  int nb, r;
  std::tie(nb, r) = computeBatchSize(n, bsize);

  // Metadata
  std::ofstream meta("smt_data/meta.txt");
  meta << 'k' << ' ' << k << std::endl;
  meta << "nb" << ' ' << nb << std::endl;
  meta.close();

  // Processing
  tbb::parallel_for(0, nb, [&](const auto i) {
    auto start = i * bsize;
    auto end = (i + 1) * bsize;
    auto *M = createDenseMT(fasta, k, start, end);
    auto *S = new arma::SpMat<uint64_t>(*M);
    delete M;

    {std::unique_lock<std::mutex> lock(mtx); save_queue.push(S);cv.notify_one();}
    
  });

  done_processing.store(true);  // Notify the saveSMT that we're done processing
  cv.notify_one();

  std::unique_lock<std::mutex> lock(mtx);
  wrt.wait(lock, [] { return done_writing.load(); });
    
  if (r != 0) {
    int start { nb * bsize };
    int end { nb * bsize + r };
    const auto *M { createDenseMT(fasta, k, start, end) };
    const auto *S { new arma::SpMat<uint64_t>(*M) };
    delete M;
    S->save(smtdb, arma::arma_binary);
    delete S;
  }

  smtdb.close();

}
