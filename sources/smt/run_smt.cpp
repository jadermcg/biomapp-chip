#include "smt.h"
#include "smt_utils.h"

#include <fstream>
#include <string>
#include <vector>
#include <lmdb.h>


int main(int argc, char* argv[]) {
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, tbb::this_task_arena::max_concurrency());
  std::string fastaPath;
  int k = 0;
  int s = 256;
  
  // Verificar se há número suficiente de argumentos
  if (argc < 7) {
    std::cerr << "Use: smt -i <fasta path> -k <size of kmer> -s <priori memory allocation>\n";
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

    else if (arg == "-s") {
      s = std::stoi(argv[i + 1]);
    }
    
    else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 1;
    }
  }
  
  // Read fasta file
  std::vector<std::string> fasta { readFasta(fastaPath) };
  std::thread process(processMT, fasta, k, s);
  std::thread save(saveMT);
  process.join();
  save.join();

  return 0;
}