#include <armadillo>
#include "smt_operations.h"
#include "smt_utils.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <fstream>

int main(int argc, char* argv[]) {
  
  int d = 0;
  std::string path2kmers = "";
  
  // Verificar se há número suficiente de argumentos
  if (argc < 5) {
    std::cerr << "Uso: hsib -kmers <kmers.txt> -d <number of mutations> < <hmap file>\n";
    return 1;
  }
  
  // Iterar através dos argumentos da linha de comando
  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    
    if (arg == "-d") {
      d = std::stoi(argv[i + 1]);
    }
    
    else if (arg == "-kmers") {
      path2kmers = argv[i + 1];
    }
    
    else {
      std::cerr << "Argumento desconhecido: " << arg << "\n";
      return 1;
    }
  }
  
  // Read map
  tbb::concurrent_hash_map <std::string, uint64_t> hmap;
  std::string kmer;
  uint64_t count;
  while (std::cin >> kmer >> count) hmap.emplace(kmer, count);
  
  // Read kmers
  std::vector<std::string> kmers = readkmers(path2kmers);
  
  // Call the hash function with the parsed arguments
  hsib(hmap, kmers, d);
  
  return 0;
}
