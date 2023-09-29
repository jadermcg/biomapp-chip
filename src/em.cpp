#include "oops.h"
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
  if (argc < 13) {
    std::cerr << "Uso: em -i <fasta> options\n   -type <oops, zoops or anr>\n   -k <size of kmer> \n   -niter <number of em iterations>\n   -cutoff <small number for convergence controll> \n   -n <number of models>\n";
    return 1;
  }
  
  std::string type = "";
  std::string path = "";
  int niter = 0;
  double cutoff = 0.0;
  int k = 0;
  int nmodels = 0;
  
  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    
    if (arg == "-type") {
      type = argv[i + 1];
    }
    
    else if (arg == "-i") {
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
    
    else if (arg == "-n") {
      nmodels = std::stoi(argv[i + 1]);
    }
    
    else {
      std::cerr << "Argumento desconhecido: " << arg << "\n";
      return 1;
    }
  }
  
  // Build siblings models
  std::vector<arma::mat> models;
  for(const auto &entry : fs::directory_iterator("smt_data/kdive_dir")) {
    std::ifstream file(entry.path());
    std::string line;
    arma::mat model(4, k);
    model.fill(1.00);
    std::string kmer;
    int count = 0;
    double total = 0.0;
    while(file >> kmer >> count) {
      int i = 0;
      for (const auto &c : kmer) {
        int line = char2int(c);
        model(line, i++) += count;
        total += count;
      }
    }
    file.close();
    models.push_back(model / (total + 4));
  }
  
  // Run EM
  const auto fasta = readFasta(path);
  arma::mat beta = createMarkovChain(fasta, 0);
  arma::mat new_alpha;
  int ret =  std::system("rm -Rf smt_data/models");
  ret = std::system("mkdir -p smt_data/models");
  
  if (type == "oops") {
    tbb::blocked_range<size_t> r(0, models.size());
    tbb::parallel_for(r, [&](const auto &r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        const auto &alpha = models[i];
        arma::mat new_alpha = oops(fasta, alpha, beta, cutoff, niter, 1.0);  // Ajuste os argumentos conforme necessário
        new_alpha.save("smt_data/models/m" + std::to_string(i + 1), arma::csv_ascii);
      }
    });
  }
  
  else if (type == "zoops") {
    tbb::blocked_range<size_t> r(0, models.size());
    tbb::parallel_for(r, [&](const auto &r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        const auto &alpha = models[i];
        arma::mat new_alpha = zoops(fasta, alpha, beta, cutoff, niter, .5);
        new_alpha.save("smt_data/models/m" + std::to_string(i+1), arma::csv_ascii);
      }
    });
  }
  
  else if (type == "anr") {
    
  }
  
  else {
    std::cerr << "Tipo inválido: precisa ser oops, zoops ou anr\n";
    return 1;
  }
  
  
  return 0;
}