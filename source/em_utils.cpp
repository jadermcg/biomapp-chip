#include "em_utils.h"
#include "utils.h"

//'Create PWM models from sibligs kmers.
//'@name build_models_from_sibligs.
//'@param k Kmer size.
//'@return A std::vector with models in arma::mat format.
std::unique_ptr<std::vector<arma::mat>> build_models_from_sibligs(const int k) {
  auto models_ptr = std::make_unique<std::vector<arma::mat>>();
  std::vector<arma::mat> &models = *models_ptr;
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
  
  return models_ptr;
}