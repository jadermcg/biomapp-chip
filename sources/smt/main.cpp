#include "smt.h"
#include "smt_utils.h"

int main(int argc, char **argv) {
    
    std::vector<std::string> fasta{ readFasta(argv[1]) };
    int len { std::stoi(argv[2]) };
    int choice { std::stoi(argv[3]) };

    int k { 30 };
    const int n = len - 0;
    const int t = fasta[0].size();
    const int m = t - k + 1;
    uint64_t nr = m * n * k;


    if (choice == 1) {
        auto *A = createDenseMT(fasta, k, 0, len);
        //arma::SpMat<uint64_t> S{*A};
        delete A;
    }

    else if (choice == 2) {
        auto *A = createDenseMT2(fasta, k, 0, len);
        delete [] A;
    }

    return 0;

}