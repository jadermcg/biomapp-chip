#include "smt_operations.h"


int main(int argc, char **argv) {

    if (argc < 3) {
        std::cerr << "Uso: ksearch -kmer <kmer>\n";
        return 1;
    }

    std::string kmer = argv[2];
    int count = ksearch(kmer);

    std::cout << kmer << ": " << count << std::endl;

    return 0;
}