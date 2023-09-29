#include "smt.h"
#include "smt_utils.h"

#include <fstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  std::string fastaPath;
  int k = 0;
  int compression = 0;
  
  // Verificar se há número suficiente de argumentos
  if (argc < 5) {
    std::cerr << "Uso: smt -i <caminho_fasta> -k <size of kmer> -c <without compression 0 or with compression 1>\n";
    return 1;
  }
  
  // Iterar através dos argumentos da linha de comando
  for (int i = 1; i < argc; i += 2) {
    std::string arg = argv[i];
    if (arg == "-i") {
      fastaPath = argv[i + 1];
    }
    
    else if (arg == "-k") {
      k = std::stoi(argv[i + 1]);
    }
    
    else if (arg == "-c") {
      compression = std::stoi(argv[i + 1]);
    }
    
    else {
      std::cerr << "Argumento desconhecido: " << arg << "\n";
      return 1;
    }
  }
  
  std::vector<std::string> fasta = readFasta(fastaPath);
  
  if (compression) createCompactMT(fasta, k); else createSparseMT(fasta, k);
  
  return 0;
}