#include "smt_operations.h"
#include "smt_utils.h"
#include <stack>
using namespace tbb;

extern std::vector<std::string> fasta;

//'search exact kmer and return your counts.
//'@name ksearch. 
//'@param kmer Kmer for search into SMT.
//'@return The number of occurrences of kmer.
int ksearch(const std::string kmer) {
  std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
  std::ifstream meta("smt_data/meta.txt");
  int k, nb;
  std::string temp;
  meta >> temp >> k;
  meta >> temp >> nb;
  meta.close();

  int count = 0;
  arma::SpMat<uint64_t> S;
  for (size_t i = 0; i < nb; ++i) {

    S.load(smtdb, arma::arma_binary);
    uint64_t node = 0;
    
    for (size_t j = 0; j < k; ++j) {
      
      int symbol = char2int(kmer[j]);
      uint64_t next = S(node, symbol);

      if (next == 0) break;
      
      node = next;
    }

    count += S(node, 4);
  }


  smtdb.close();
  return count;
}

//'Computes hamming distance efficiently.
//'@name hDist.
//'@param str1 First string to compare.
//'@param str2 Second string to compare.
//'@return Hamming distance between str1 and str2.
int hDist(const std::string &str1, const std::string &str2) {
  return std::inner_product(str1.begin(), str1.end(), str2.begin(), 0, std::plus<>(), std::not_equal_to<>());
}

//'Kmer direct Search.
//'@name busca_direta.
//'@param fasta Dataset of sequences.
//'@param kmers kmers to find sibligs.
//'@param d Number of mutations allowed.
//'@return C++ String HashMap of siblings of kmers.
std::map<std::string, std::map<std::string, int>> busca_direta(std::vector<std::string> &fasta, std::vector<std::string> &kmers, int d) {
  int n = fasta.size();
  int t = fasta[0].size();
  int k = kmers[0].size();
  int m = t - k + 1;
  std::map<std::string, std::map<std::string, int>> hmap;

  for (int i = 0; i < n; ++i) {
    std::string seq = fasta[i];
    for (int j = 0; j < m; ++j) {
      std::string km = seq.substr(j, k);
      for (const auto &kmer : kmers) {
        int e = hDist(kmer, km);
        if (e <= d) {
          hmap[kmer][km] += 1;
        }
      }
    }
  }

  return hmap;
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
  std::ifstream file("smt_data/meta.txt");
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

//'Auxiliary Recursive kdive function.
//'@name kdive_.
//'@param M SMT data.
//'@param hmap C++ String HashMap for store sibligs.
//'@param kmer Kmer for search siblings.
//'@param k Size of kmer.
//'@param d Number of mutations allowed.
//'@param node Root of SMT.
//'@param l Current number of mutations. Needs to be less than k.
//'@param j Current index of kmer.
void kdive_(const arma::Mat<uint64_t> &M, tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> &hmap, const std::string &kmer, const int &k, const int &d, int node, int l, int j) {
  
  if (j == k) {
    tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>>::accessor outer_acc;
    tbb::concurrent_hash_map<std::string,uint64_t>::accessor inner_acc;

    if (hmap.insert(outer_acc, kmer)) outer_acc->second = tbb::concurrent_hash_map<std::string, uint64_t>();
    
    if (outer_acc->second.insert(inner_acc, index2kmer(M(node, 5), k))) inner_acc->second = 0;
    inner_acc->second += M(node, 4);
    
    inner_acc.release();
    outer_acc.release();
    return;
  }
  
  tbb::parallel_for(0, 4, 1, [&](size_t i) {
    int next = M(node, i);
    if (next != 0) {
      char c = int2char(i);
      int hd = (kmer[j] == c) ? 0 : 1;
      if (l + hd <= d) {
        kdive_(M, hmap, kmer, k, d, next, l + hd, j + 1);
      }
    }
  });
}

//'Auxiliary Iterative kdive function.
//'@name kdive_.
//'@param M SMT data.
//'@param hmap C++ String HashMap for store sibligs.
//'@param kmer Kmer for search siblings.
//'@param k Size of kmer.
//'@param d Number of mutations allowed.
void kdive_iterativo(const arma::Mat<uint64_t> &M, tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> &hmap, const std::string &kmer, const int &k, const int &d) {
    struct Frame {
        int node, l, j;
    };
    
    std::stack<Frame> pilha;
    pilha.push({0, 0, 0});  // Inicializa a pilha com o estado inicial

    while (!pilha.empty()) {
        Frame frameAtual = pilha.top();
        pilha.pop();

        int node = frameAtual.node;
        int l = frameAtual.l;
        int j = frameAtual.j;

        if (j == k) {
            tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>>::accessor outer_acc;
            tbb::concurrent_hash_map<std::string,uint64_t>::accessor inner_acc;

            if (hmap.insert(outer_acc, kmer)) outer_acc->second = tbb::concurrent_hash_map<std::string, uint64_t>();
            
            if (outer_acc->second.insert(inner_acc, index2kmer(M(node, 5), k))) inner_acc->second = 0;
            inner_acc->second += M(node, 4);
            
            inner_acc.release();
            outer_acc.release();
        } else {
            for (size_t i = 0; i < 4; ++i) {
                int next = M(node, i);
                if (next != 0) {
                    char c = int2char(i);
                    int hd = (kmer[j] == c) ? 0 : 1;
                    if (l + hd <= d) {
                        pilha.push({next, l + hd, j + 1});  // Empilha o próximo estado
                    }
                }
            }
        }
    }
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
  std::ifstream meta("smt_data/meta.txt");
  std::string temp;
  int k, nb;
  meta >> temp >> k; // read first line
  meta >> temp >> nb; // read second line
  meta.close();
  
  std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
  arma::SpMat<uint64_t> S;
  arma::Mat<uint64_t> M;
  for (size_t i = 0; i < nb; ++i) {
    S.load(smtdb);

    M = S;
    for (const auto &kmer : kmers) kdive_(M, hmap, kmer, k, d,0,0,0);
  }
  
  smtdb.close();
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
      const auto &km = it->first;
      const auto &count = it->second;
      int hd = hDist(km, kmer);
      if (hd <= d) {
        fhsib << km << ' ' << count << std::endl;
      }
    }
    
    fhsib.close();
    
  });
}




