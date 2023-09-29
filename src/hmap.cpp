#include <armadillo>
#include "smt_operations.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <fstream>

int main(int argc, char* argv[]) {
  
  // Verificar se há número suficiente de argumentos
  if (argc < 1) {
    std::cerr << "Uso: hmap > <output file.txt>\n";
    return 1;
  }
  
  // Call the hash function with the parsed arguments
  const auto &hash = hmap();
  
  // Imprimir o vector ordenado
  for (const auto &pair : hash) {
    std::cout << pair.first << " " << pair.second << '\n';
  }
  
  return 0;
}
