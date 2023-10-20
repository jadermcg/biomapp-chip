#include <armadillo>
#include "smt_operations.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <fstream>
#include <tbb/concurrent_hash_map.h>

int main(int argc, char* argv[]) {
  
  int k = 0;
  
  // Verificar se há número suficiente de argumentos
  if (argc < 3) {
    std::cerr << "Uso: khmap -k <size of kmer>\n";
    return 1;
  }
  
  // Iterar através dos argumentos da linha de comando
  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    
    
    if (arg == "-k") {
      k = std::stoi(argv[i + 1]);
    }
    
    else {
      std::cerr << "Invalid argument: " << arg << "\n";
      return 1;
    }
  }
  
  // Chamar a função de hash com os argumentos analisados
  tbb::concurrent_hash_map<std::string, uint64_t> hmap { khmap(k) };
  
  // Imprimir o vetor ordenado
  for (const auto &pair : hmap) {
    std::cout << pair.first << " " << pair.second << '\n';
  }
  
  return 0;
}