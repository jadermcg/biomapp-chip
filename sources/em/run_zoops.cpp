#include "zoops.h"
#include "utils.h"
#include "prob_utils.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  
  // Verificar se há número suficiente de argumentos
  if (argc < 9) {
    std::cerr << "Uso: oops -i <fasta> options\nOptions:\n   -k <size of kmer>\n   -niter <number of em iterations>\n   -cutoff <small number for convergence controll>\n";
    return 1;
  }
  
  std::string path = "";
  int niter = 0;
  double cutoff = 0.0;
  int k = 0;
  
  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    
    if (arg == "-i") {
      path = argv[i + 1];
    }
    
    else if (arg == "-niter") {
      niter = std::stoi(argv[i + 1]);
    }
    
    else if (arg == "-cutoff") {
      cutoff = std::stod(argv[i + 1]);
    }
    
    else if (arg == "-k") {
      k = std::stoi(argv[i + 1]);
    }
    
    else {
      std::cerr << "Argumento desconhecido: " << arg << "\n";
      return 1;
    }
  }
  
  // Run oops EM
  const auto fasta = readFasta(path);
  arma::mat beta = createMarkovChain(fasta, 0);
  arma::mat alpha = fasta2alpha(fasta, k);
  int ret =  std::system("rm -Rf zoops/models");
  ret = std::system("mkdir -p zoops/models");
  arma::mat new_alpha = zoops(fasta, alpha, beta, cutoff, niter, 1.0);  // Ajuste os argumentos conforme necessário
  new_alpha.save("zoops/models/m1", arma::csv_ascii);
  
  return 0;
}