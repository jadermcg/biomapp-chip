#include "smt.h"
#include "smt_utils.h"

//'@keywords internal
std::tuple<int, int, int> computeBatchSize(int n, int k, int t) {
  int bsize = 0;
  int nb = 0;
  int r = 0;
  
  if (n * t > 1e6) {
    bsize = 1e6 / t / std::log(k);
    nb = n / bsize;
    r = n - nb * bsize;
  }
  else {
    nb = 1;
    bsize = n;
  }
  
  return std::make_tuple(bsize, nb, r);
}

//'This function shows the sparsity level of the SMT matrix.
//'@name sparsity_level
//'@param SMT The sparse Motif Tree.
//'@return The sparsity level of SMT.
//[[Rcpp::export]]
double sparsity_level() {
  // double num_non_zeros = static_cast<double>(SMT.n_nonzero);
  // double total_elements = static_cast<double>(SMT.n_rows * SMT.n_cols);
  // return 1.0 - (num_non_zeros / total_elements);
  return -1;
}

//'@keywords internal
arma::Mat<uint64_t>* createDenseMT(const std::vector<std::string> &fasta, const int k, const int start, const int end) {
 const int n = end - start;
 const int t = fasta[0].size();
 const int m = t - k + 1;
 
 arma::uword nr = m * n * k;
 auto *MT = new arma::Mat<uint64_t>(nr, 6);
 int next = 0;
 
 for (int i = start; i < end; ++i) {
   const auto &seq = fasta[i];
   
   for (int j = 0; j < m; ++j) {
     uint64_t node = 0;
     uint64_t index = 0;
     
     for (int l = 0; l < k; ++l) {
       char c = seq[j + l];
       int symbol = char2int(c);
       index = index * 4 + symbol;
       int current = MT->at(node, symbol);
       
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

//'Creates SMT matrix from the sequences.
//'Creates SMT matrix from the sequences.
//'@name createSMT
//'@param fasta The Dataset of sequences.
//'@param k The Size of kmers.
//'@return SMT The Sparse Motif Tree.
//[[Rcpp::export]]
void createSparseMT(const std::vector<std::string> &fasta, const int k) {
 
 int n = fasta.size();
 int t = fasta[0].size();
 
 // Number of batches
 int bsize, nb, r;
 std::tie(bsize, nb, r) = computeBatchSize(n, k, t);
 
 //setupBuffer
 int ret = std::system("rm -Rf smt_data");
 ret = std::system("mkdir smt_data");
 
 // Saving metadata
 std::ofstream meta("smt_data/metadata.txt");
 meta << "k " << k << std::endl;
 meta << "nb " << ((r == 0) ? nb : nb + 1) << std::endl;
 meta.close();
 
 // Fluxo de saída
 std::ofstream smtdb("smt_data/SMT.db", std::ios::binary);
 std::mutex mtx;
 tbb::parallel_for(0, nb, 1, [&](size_t i) {
   
   int start = i * bsize;
   int end = (i + 1) * bsize;
   const auto *M = createDenseMT(fasta, k, start, end);
   const auto *S = new arma::SpMat<uint64_t>(*M);
   
   // Critical
   {std::lock_guard<std::mutex> lock(mtx); S->save(smtdb, arma::arma_binary);}
   
   delete M;
   delete S;
   
 });
 
 // IF threre are remainder
 if (r != 0) {
   int start = nb * bsize;
   int end =  nb * bsize + r;
   const auto *M = createDenseMT(fasta, k, start, end);
   const auto *S = new arma::SpMat<uint64_t>(*M);
   S->save(smtdb, arma::arma_binary);
   delete M;
   delete S;
 }
 
 smtdb.close();
}

//'Creates SMT matrix from the sequences.
//'Creates SMT matrix from the sequences.
//'@name createSMT
//'@param fasta The Dataset of sequences.
//'@param k The Size of kmers.
//'@return SMT The Sparse Motif Tree.
//[[Rcpp::export]]
void createCompactMT(const std::vector<std::string> &fasta, const int k) {
  
  int n = fasta.size();
  int t = fasta[0].size();
  int m = t - k + 1;
  
  // Number of batches
  int bsize, nb, r;
  std::tie(bsize, nb, r) = computeBatchSize(n, k, t);
  
  //setupBuffer
  int ret = std::system("rm -Rf smt_data");
  ret = std::system("mkdir smt_data");
  
  // Saving metadata
  std::ofstream meta("smt_data/metadata.txt");
  meta << "k " << k << std::endl;
  meta << "nb " << ((r == 0) ? nb : nb + 1) << std::endl;
  
  // Fluxo de saída
  std::ofstream smtdb("smt_data/SMT.db", std::ios::binary);
  std::stringstream ss;
  int maxDstSize = LZ4_compressBound(m*bsize*k*6);
  std::vector<char> compressed(maxDstSize);
  std::mutex mtx;
  tbb::parallel_for(0, nb, 1, [&](size_t i) {
    
    int start = i * bsize;
    int end = (i + 1) * bsize;
    const auto *M = createDenseMT(fasta, k, start, end);
    
    // Critical
    {
      std::lock_guard<std::mutex> lock(mtx);
      M->save(ss, arma::arma_binary);
      std::string uncompressedData = ss.str();
      int compressedSize = LZ4_compress_fast(uncompressedData.data(), compressed.data(), uncompressedData.size(), maxDstSize, 1);
      meta << compressedSize << " " << uncompressedData.size() << std::endl;
      smtdb.write(compressed.data(), compressedSize);
      ss.str("");
      ss.clear();
    }
    
    delete M;
  });
  
  // IF threre are remainder
  if (r != 0) {
    int start = nb * bsize;
    int end =  nb * bsize + r;
    const auto *M = createDenseMT(fasta, k, start, end);
    M->save(ss, arma::arma_binary);
    std::string uncompressedData = ss.str();
    int compressedSize = LZ4_compress_fast(uncompressedData.data(), compressed.data(), uncompressedData.size(), maxDstSize, 1);
    meta << compressedSize << " " << uncompressedData.size() << std::endl;
    smtdb.write(compressed.data(), compressedSize);
    ss.str("");
    ss.clear();
    delete M;
  }
  
  meta.close();
  smtdb.close();
  
}
