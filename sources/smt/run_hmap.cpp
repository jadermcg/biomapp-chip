#include "hmap.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <sqlite3.h>

extern tbb::concurrent_unordered_map<std::string, uint64_t> final_hash;

int main(int argc, char* argv[]) {

  int top = 0;

  // Verificar se há número suficiente de argumentos
  if (argc < 1) {
    std::cerr << "Uso: hmap\n";
    return 1;
  }


  // Call the hash function with the parsed arguments
  std::ifstream meta("smt_data/meta.txt");
  int nb, k;
  std::string str;
  meta >> str >> k;
  meta >> str >> nb;
  meta.close();

  hmap(nb, k);
  std::thread c1(consumer);
  std::thread c2(consumer);
  std::thread c3(consumer);

  c1.join();
  c2.join();
  c3.join();
  
  
  // Saving hash in file
  std::mutex mtx;
  std::ofstream outfile("smt_data/hmap.txt");
  tbb::parallel_for(final_hash.range(), [&](auto &range)  {
    std::ostringstream buffer;
    for (auto it { range.begin() }; it != range.end(); ++it) {
      buffer << it->first << " " << it->second << std::endl;

      if (buffer.tellp() >= 16384) {
        std::unique_lock<std::mutex> lock(mtx);
        outfile << buffer.str();
        buffer.str("");
        buffer.clear();
      }
    }

    if (!buffer.str().empty()) {
      std::unique_lock<std::mutex> lock(mtx);
      outfile << buffer.str();
    }
  });

  outfile.close();

  return 0;
}
