#include "smt_operations.h"
#include "smt_utils.h"

int main(int argc, char **argv) {

    if (argc < 7) {
        std::cerr << "Use: dsearch -i fasta <path to fasta> -kmers <path to kmers> -d <number of mutations>\n";
        return -1;
    }

    std::string path;
    std::string path2kmers;
    int d;

    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];

        if (arg == "-i") {
            path = argv[i + 1];
        }

        else if (arg == "-kmers") {
            path2kmers = argv[i + 1];
        }

        else if (arg == "-d") {
            d = std::stoi(argv[i + 1]);
        }

        else {
            std::cerr << "Invalid parameter"<< std::endl;
            return -1;
        }
    }

    std::vector<std::string> fasta = readFasta(path);
    std::vector<std::string> kmers = readkmers(path2kmers);

    std::map<std::string, std::map<std::string, int>> hmap = busca_direta(fasta, kmers, d);

    int ret = std::system("rm -Rf smt_data/dsearch_dir");
    ret = std::system("mkdir -p smt_data/dsearch_dir");
    for (const auto &kmer : hmap) {
        std::ofstream f("smt_data/dsearch_dir/" + kmer.first + ".txt");    
        for (const auto &km : kmer.second) {
            f << km.first << " " << km.second << std::endl;
        }
        f.close();
    }

    return 0;
}