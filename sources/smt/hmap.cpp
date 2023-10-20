#include "hmap.h"
#include "smt_utils.h"

extern std::vector<std::string> fasta;
std::ifstream smtdb("smt_data/SMT.db", std::ios::binary);
tbb::concurrent_unordered_map<std::string, uint64_t> final_hash;
std::mutex mtx;
std::atomic_bool still_processing(true);
tbb::concurrent_queue<tbb::concurrent_unordered_map<std::string, uint64_t>> maps;
#define mem_reserv 10000

void consumer() {
  while (!maps.empty() || still_processing.load()) {
      tbb::concurrent_unordered_map<std::string, uint64_t> map;
      map.reserve(mem_reserv);
      while (maps.try_pop(map)) {
          for(auto& kv : map) {
              final_hash[kv.first] += kv.second;
          }
      }
  }
}

//'Compute fast hashmap from smt_data.
//'@name fast_hash
//'@param path Path to smt_data.
//'@param nthreads Number of threads.
//'@return C++ String Hash Map.
void hmap(int nb, int k) {

  tbb::parallel_for(tbb::blocked_range<size_t>(0, nb), [&] (const auto &r) {
    tbb::concurrent_unordered_map<std::string, uint64_t> hash;
    hash.reserve(mem_reserv);

    for (size_t i = r.begin(); i < r.end(); ++i) {

      // Load file
      auto *S { new arma::SpMat<uint64_t>() };
      {std::unique_lock lock(mtx); S->load(smtdb, arma::arma_binary);}

      // Processing
      arma::Col<uint64_t> counts{ S->col(4) };
      arma::Col<uint64_t> kmers { S->col(5) };
      delete S;
      arma::uvec nonZeroIndices { arma::find(counts > 0) };
      counts = counts(nonZeroIndices);
      kmers = kmers(nonZeroIndices);

      for (size_t i {0}; i < nonZeroIndices.n_elem; ++i) {
        std::string kmer { index2kmer(kmers[i], k) };
        hash[kmer] += counts[i];
      }
    }
    maps.push(hash);
  });

  still_processing.store(false);

}