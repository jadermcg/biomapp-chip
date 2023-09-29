#include <armadillo>
#include "smt_operations.h"
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
    std::cerr << "Uso: kdive -kmers <path to kmers> -d <number of mutations>\n";
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
  
  // Read kmers
  std::vector<std::string> kmers;
  std::ifstream file(path2kmers);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      kmers.push_back(line);
    }
    file.close();
  } else {
    std::cerr << "Não foi possível abrir o arquivo!" << std::endl;
  }
  
  // Call the hash function with the parsed arguments
  const auto &hmap = kdive(kmers, d);
  
  // Criando diretório de saída
  int ret = system("rm -Rf smt_data/kdive_dir");
  mkdir("smt_data/kdive_dir", 0777);
  
  // Iterando pelo concurrent_hash_map externo
  //for (auto it = hmap.begin(); it != hmap.end(); ++it) {
  tbb::parallel_for_each(hmap.begin(), hmap.end(), [&](const auto& it) {
    const auto &kmer = it.first;
    const auto &map = it.second;
    std::ofstream fdive("smt_data/kdive_dir/" + kmer + ".txt");
    
    // Iterando pelo concurrent_hash_map interno
    for (auto inner_it = map.begin(); inner_it != map.end(); ++inner_it) {
    //tbb::parallel_for_each(map.begin(), map.end(), [&](const auto& inner_it) {
      const auto &sibling = inner_it->first;
      const auto &counts = inner_it->second;
      fdive << sibling << ' ' << counts << std::endl;
    }
    
    fdive.close();
  
  });
  
  
  return 0;
}
