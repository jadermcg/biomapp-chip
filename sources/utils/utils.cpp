#include "utils.h"

//'Read all file names into a dirPath.
//'getFilenames
//'@param dirPath Directory path.
//'@return C++ string vector with directory file names.
//[[Rcpp::export]]
std::vector<std::string> getFilenames(const std::string& dirPath) {
 std::vector<std::string> filenames;
 boost::filesystem::path path(dirPath);
 
 try {
   if (boost::filesystem::exists(path) && boost::filesystem::is_directory(path)) {
     for (boost::filesystem::directory_iterator itr(path); itr != boost::filesystem::directory_iterator(); ++itr) {
       filenames.push_back(itr->path().filename().string());
     }
   } else {
     std::cerr << "Caminho não encontrado ou não é um diretório\n";
   }
 } catch (const boost::filesystem::filesystem_error& e) {
   std::cerr << "Erro ao acessar o diretório: " << e.what() << '\n';
 }
 
 return filenames;
}

//'Read fasta dataset.
//'@name readFasta
//'@param filepath Path to fasta dataset.
//'@return C++ vector string whth each line is a sequence in fasta dataset.
std::vector<std::string> readFasta(const std::string& filepath) {
  std::ifstream file(filepath);
  std::vector<std::string> data;
  std::string line;
  std::string seq;
  
  int i = 0;
  while (std::getline(file, line)) {
    
    if (line[0] == '>' && i != 0) {
      if (seq.find('N') == std::string::npos) data.push_back(seq);
      seq.clear();
    }
    
    else if (line[0] != '>' && i != 0) {
      seq += line;
    }
    
    ++i;
  }
  
  if (seq.find('N') == std::string::npos) data.push_back(seq);
  
  return data;
}

//'Converts char nucleotide A,C,G,T in int 0,1,2,3.
//'@name char2int
//'@param c char to convert for.
//'@return The correspondent int or throws a invalid character.
int char2int(char c) {
   switch(c) {
     case 'A': return 0;
     case 'C': return 1;
     case 'G': return 2;
     case 'T': return 3;
     case 'X': return -1;
     default: throw std::invalid_argument("Invalid character");
   }
 }

//'Converts int 0,1,2,3 to nucleotide char A,C,G,T.
//'@name int2char
//'@param i int to convert for.
//'@return The correpondent nucleotide character or throws a invalid integer.
char int2char(int i) {
   switch(i) {
     case 0: return 'A';
     case 1: return 'C';
     case 2: return 'G';
     case 3: return 'T';
     case -1: return 'X';
     default: throw std::invalid_argument("Invalid integer");
   }
 }

//'Convert a kmer into a corresponding integer.
//'@name kmer2index
//'@param kmer Kmer to convert for.
//'@return A corresponding index to kmer.
int kmer2index(const std::string &kmer) {
   int k = kmer.size();
   int index = 0;
   for (int i = 0; i < k; ++i) {
     index = index * 4 + char2int(kmer[i]);
   }
   return index;
 }

//'Convert a integer to corresponding kmer.
//'@name index2kmer
//'@param index Index to converting for.
//'@param k The size of kmer.
//'@return The corresponding kmer.
std::string index2kmer(int index, int k) {
   std::string kmer(k, 'A'); // Inicializa a string kmer com tamanho k e todos os caracteres como 'A'
   
   for (int i = k - 1; i >= 0; --i) {
     int residue = index % 4;
     kmer[i] = int2char(residue);
     index /= 4;
   }
   
   return kmer;
 }

//'Computes fast log2 function.
//'@name fastlog2
//'@param x A positive double value.
//'@return Fast log2 of x.
double fastlog2(double x) {
   union {
   double d;
   int64_t i;
 } conv;
   
   conv.d = x;
   int64_t exponent = ((conv.i >> 52) & 0x7FF) - 1023;
   conv.i = (conv.i & ~(static_cast<int64_t>(0x7FF) << 52)) | (static_cast<int64_t>(0x3FF) << 52);
   return (exponent + conv.d) - 1.0;
 }

//'Computes a log fast log function.
//'@name fastlog.
//'@param x A positive double value.
//'@return Fast Log of x.
double fastlog(double x) {
   return 0.69314718f * fastlog2(x);
 }

//'Computes hamming distance for two strings.
//'@name hamming_distance_parallel
//'@param str1 The first string.
//'@param str2 The second string.
//'@return The hamming distance between str1 and str2.
int hamming_distance_parallel(const std::string &str1, const std::string &str2) {
   if (str1.size() != str2.size()) {
     throw std::invalid_argument("Strings must have the same length to compute Hamming distance.");
   }
   
   int distance = 0;
   size_t size = str1.size();
   #pragma omp parallel for reduction(+:distance)
   for (size_t i = 0; i < size; ++i) {
     if (str1[i] != str2[i]) {
       ++distance;
     }
   }
   
   return distance;
 }

//'This function computes the overlap string between two string a and b.
//'@name string_corr
//'@param a First string.
//'@param b Second string.
//'@return The correlation string between both strings a and b.
//[[Rcpp::export]]
std::string corr(const std::string &a, const std::string b) {
   int k = a.size();
   std::string c;
   for (int i = 0; i < k; ++i) {
     bool isover = true;
     for (int j = i; j < k; ++j) {
       if (a[j] != b[j-i]) {
         isover = false;
         break;
       }
     }
     if (isover) c[i] = '1';
     else c[i] = '0';
   }
   return c;
 }

//'Computes the correlation frequency from correlation string.
//'@name corr_freq
//'@param correlation_str Correlation string computed from corr.
//'@return Frenquency correlation.
//[[Rcpp::export]]
double corr_freq(const std::string &correlation_str) {
   double corr_f = 0.0;
   int k = correlation_str.size();
   
   for (int i = 0; i < k; ++i) {
     char c = correlation_str[i];
     double x = std::stod(std::string(1, c));
     corr_f += x * pow(0.5, i);
   }
   
   return corr_f;
 }

//'Computes the correlation frequency between two strings.
//'@name fast_corr_freq
//'@param a First string.
//'@param b Second string.
//'@return Frenquency correlation.
//[[Rcpp::export]]
double fast_corr_freq(const std::string &a, const std::string b) {
     int k = a.size();
     double corr_f = 0.0;
     
     for (int i = 0; i < k; ++i) {
       bool isover = true;
       for (int j = i; j < k; ++j) {
         if (a[j] != b[j-i]) {
           isover = false;
           break;
         }
       }
       if (isover) corr_f += pow(0.5, i);
     }
     return corr_f;
   }

//'Compute DKL between alpha and beta models. 
//'@name computeDKL
//'@param alpha PWM model.
//'@param beta 0-order Markov model.
//'@param kmer Target sequence.
//'@return Kullback-Leibler divergence between alpha and beta models.
//[[Rcpp::export]]
double computeDKL(const arma::mat &alpha, const arma::mat &beta, const std::string &kmer) {
  int k = alpha.n_cols;
  double dkl = 0.0;
  
  for (int i = 0; i < k; ++i) {
    int idx = char2int(kmer[i]);
    dkl += alpha(idx,i) * (std::log(alpha(idx,i)) - std::log(beta(idx)));
  }
  
  return dkl;
}

//'Compute DKL between alpha and uniform models. 
//'@name computeDKLU
//'@param alpha PWM model.
//'@param kmer Target sequence.
//'@return Kullback-Leibler divergence between alpha and uniform models.
//[[Rcpp::export]]
double computeDKLU(const arma::mat &alpha, const std::string &kmer) {
  int k = alpha.n_cols;
  double dkl = 0.0;
  
  for (int i = 0; i < k; ++i) {
    int idx = char2int(kmer[i]);
    dkl += alpha(idx,i) * (2 + std::log(alpha(idx,i)));
  }
  
  return dkl;
}

//'Compute score from PWM model. 
//'@name score
//'@param alpha PWM model.
//'@return Simple score from PWM model.
//[[Rcpp::export]]
int computeSCORE(const arma::mat &alpha) {
  arma::rowvec max_values = arma::max(alpha, 0);
  int score = arma::sum(1 - max_values);
  return score;
}
