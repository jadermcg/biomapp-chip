#include "smt_operations.h"
#include "smt_utils.h"

using namespace tbb;

//'Computes hamming distance efficiently.
//'@name hDist.
//'@param str1 First string to compare.
//'@param str2 Second string to compare.
//'@return Hamming distance between str1 and str2.
int hDist(const std::string &str1, const std::string &str2) {
  return std::inner_product(str1.begin(), str1.end(), str2.begin(), 0, std::plus<>(), std::not_equal_to<>());
}

//'Compute fast hashmap from smt_data.
//'@name fast_hash
//'@param path Path to smt_data.
//'@param nthreads Number of threads.
//'@return C++ String Hash Map.
tbb::concurrent_hash_map<std::string, uint64_t> hmap() {
  tbb::concurrent_hash_map<std::string, uint64_t> hash;
  
  // Read metadata
  std::ifstream meta("smt_data/metadata.txt");
  int k, nb;
  std::string str;
  meta >> str >> k;
  meta >> str >> nb;
  meta.close();
  
  // Read smt data
  std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
  std::mutex mtx;
  tbb::parallel_for(0, nb , 1, [&](size_t i) {
    arma::SpMat<uint64_t> M;
    
    // Critical
    {std::lock_guard lock(mtx); M.load(smtdb);}
    
    arma::Col<uint64_t> counts(M.col(4));
    arma::Col<uint64_t> kmers(M.col(5));
    arma::uvec nonZeroIndices = arma::find(counts > 0);
    counts = counts(nonZeroIndices);
    kmers = kmers(nonZeroIndices);
    for (size_t i = 0; i < nonZeroIndices.n_elem; ++i) {
      std::string kmer = index2kmer(kmers[i], k);
      tbb::concurrent_hash_map<std::string, uint64_t>::accessor acc;
      hash.insert(acc, kmer);
      acc->second += counts[i];
    }
  });
  
  smtdb.close();
  
  return hash;
}

//'Count kmers in a SMT tree.
//'@name count_kmers.
//'@param S SMT data.
//'@param hmap C++ String Hash Map.
//'@param kmer Kmer for counting.
//'@param kmax Size of kmer.
//'@param j Start index of kmer. Does not need be 0.
//'@param node Current node. Does not be root.
//'@return hmap[kmer] += count.
void count_kmers(const arma::SpMat<uint64_t> &S, concurrent_hash_map<std::string, uint64_t> &hmap, const std::string &kmer, const int kmax, int j, int node) {
  
  if (j == kmax) {
    int count = S(node, 4);
    
    // Atualização thread-safe usando TBB
    concurrent_hash_map<std::string, uint64_t>::accessor acc;
    hmap.insert(acc, kmer);
    acc->second += count;
    
    return;
  }
  
  // Executa as iterações em paralelo usando TBB
  parallel_for(0, 4, 1, [&](size_t i) {
    int next = S(node, i);
    if (next > 0) {
      count_kmers(S, hmap, kmer, kmax, j + 1, next);
    }
  });
}

//'Auxiliary function to hash.
//'@name hash_.
//'@param S SMT data.
//'@param hmap C++ String Hash Map.
//'@param kmer Kmer for hashing.
//'@param kmax kmax Max size of kmer.
//'@param k Size of the kmer of interest.
//'@oaram j Start index of kmer. Does not need be 0.
//'@param node Current node. Does not be root.
//'@nthreads Number os threads.
//'@return Calls function count_kmers.
void hash(const arma::SpMat<uint64_t> &S, concurrent_hash_map<std::string, uint64_t> &hmap, std::string kmer, const int kmax, const int k, int j, int node) {
  
  if (j == k) {
    count_kmers(S, hmap, kmer, kmax, j, node);
    return;
  }
  
  // Executa as iterações em paralelo usando TBB
  parallel_for(0, 4, 1, [&](size_t i) {
    int next = S(node, i);
    if (next > 0) {
      hash(S, hmap, kmer + int2char(i), kmax, k, j + 1, next);
    }
  });
}

//'Hash and count kmers.
//'@name khmap.
//'@param k Size of the kmer of interest.
//'@param path Path to SMT data.
//'@nthreads Number os threads.
//'@return C++ String HashMap of kmers and your counts.
tbb::concurrent_hash_map<std::string, uint64_t> khmap(const int k) {
  concurrent_hash_map<std::string, uint64_t> hmap;
  
  // Ler metadados
  int kmax, nb;
  std::string temp;
  std::ifstream file("smt_data/metadata.txt");
  file >> temp >> kmax;
  file >> temp >> nb;
  file.close();
  
  if (k > kmax) {
    throw std::runtime_error("K precisa ser menor que kmax!");
  }
  
  
  // Executa em paralelo usando TBB
  std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
  std::mutex mtx;
  parallel_for(0, nb, 1, [&](size_t i) {
    arma::SpMat<uint64_t> S;
    
    {std::lock_guard lock(mtx); S.load(smtdb, arma::arma_binary);}
    
    
    hash(S, hmap, "", kmax, k, 0, 0);
  });
  
  smtdb.close();
  
  return hmap;
}

//'Auxiliary kdive function.
//'@name kdive_.
//'@param S SMT data.
//'@param hmap C++ String HashMap for store sibligs.
//'@param kmer Kmer for search siblings.
//'@param k Size of kmer.
//'@param d Number of mutations allowed.
//'@param node Root of SMT.
//'@param l Current number of mutations. Needs to be less than k.
//'@param j Current index of kmer.
//'@param sibling Sibling of a kmer. 
void kdive_(const arma::SpMat<uint64_t> &S, tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> &hmap, const std::string &kmer, const int k, const int d, int node, int l, int j, std::string sibling) {
  
  if (j == k) {
    tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>>::accessor outer_acc;
    tbb::concurrent_hash_map<std::string,uint64_t>::accessor inner_acc;
    if (hmap.insert(outer_acc, kmer)) outer_acc->second = tbb::concurrent_hash_map<std::string, uint64_t>();
    if (outer_acc->second.insert(inner_acc, sibling)) inner_acc->second = 0;
    inner_acc->second += S(node, 4);
    inner_acc.release();
    outer_acc.release();
    return;
  }
  
  //for (int i = 0; i < 4; ++i) {
  tbb::parallel_for(0, 4, 1, [&](size_t i) {
    int next = S(node, i);
    if (next != 0) {
      char c = int2char(i);
      int hd = (kmer[j] == int2char(i)) ? 0 : 1;
      if (l + hd <= d) {
        kdive_(S, hmap, kmer, k, d, next, l + hd, j + 1, sibling + c);
      }
    }
  });
}

//'Search all siblings of a kmers. Do not use this function! Use hsib or fast_hsib.
//'@name kdive.
//'@param kmers List of kmers for search siblings.
//'@param d Number of mutations allowed.
//'@param path Path to SMT data.
//'@param nthreads Numbers of threads.
//'@return C++ String HashMap of siblings of kmers.
tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> kdive(const std::vector<std::string> &kmers, const int d) {
  
  tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> hmap;
  
  // Open and read metadata
  std::ifstream meta("smt_data/metadata.txt");
  std::string temp;
  int k, nb;
  meta >> temp >> k; // read first line
  meta >> temp >> nb; // read second line
  meta.close();
  
  std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
  arma::SpMat<uint64_t> S;
  for (size_t i = 0; i < nb; ++i) {
    S.load(smtdb);
    for (const auto &kmer : kmers) kdive_(S, hmap, kmer, k, d, 0, 0, 0, "");
  }
  
  return hmap;
}

//'Search siblings of kmers in a hmap data.
//'@name fast_hsib.
//'@param hmap HashMap of kmers and yours counts.
//'@param kmers Kmers to search siblings.
//'@param d Number of mutations allowed.
//'@param nthread Number os threads.
//'@return C++ String HashMap of siblings of kmers.
void hsib(const tbb::concurrent_hash_map <std::string, uint64_t> &hmap, const std::vector<std::string> &kmers, const int d) {
  
  int ret = system("rm -Rf smt_data/hsib_dir");
  mkdir("smt_data/hsib_dir", 0777);
  tbb::parallel_for(size_t(0), kmers.size(), [&](size_t i) {
    std::string kmer = kmers[i];
    std::ofstream fhsib("smt_data/hsib_dir/" + kmer + ".txt");
    
    for (auto it = hmap.begin(); it != hmap.end(); ++it) {
    //tbb::parallel_for_each(hmap.begin(), hmap.end(), [&](const auto& it) {
      const auto &km =it->first;
      const auto &count = it->second;
      int hd = hDist(km, kmer);
      if (hd <= d) {
        fhsib << km << ' ' << count << std::endl;
      }
    }
    
    fhsib.close();
    
  });
}




