#include "prob_utils.h"
#include "utils.h"
std::mt19937 gen(std::random_device{}()); 

//'Check if EM is converged.
//'@name hasConverged
//'@param cutoff If change is less than cutoff, EM has converged.
//'@param niter
//'@param alpha
//'@param convergence
//'@param changes
//'@return True or False dependent on whether the EM has converged or not.
bool hasConverged(const double cutoff, const int niter, const arma::mat &alpha, std::vector<double> &convergence, std::vector<double> &changes) {
  convergence.push_back(computeICU(alpha));
  double change = step_size(convergence);
  //std::cout << change << std::endl;
  changes.push_back(change);
  if (niter == 0 || change < cutoff) return true;
  
  return false;
}

//'Compute the step-size for EM convergence.
//'@name step_size
//'@param convergence Convergence vector.
//'@return The abs of difference between the 2 last positions of convercence vector.
double step_size(std::vector<double> convergence) {
   int n = convergence.size();
   double change = std::abs(convergence[n - 2] - convergence[n - 1]);
   return change;
}

//'Computes de complete Log-Likelihood of the data and parameters.
//'@name Q
//'@param fasta Dataset of sequences.
//'@param alpha PWM motif model.
//'@param beta Markov backround model.
//'@return The Log-Likelihood of all sequences and parameters.
double Q(const std::string &mod, const std::vector<std::string> &fasta, const arma::mat &alpha, const arma::mat &beta, const double w) {
   double q = 0.0;
   int k = alpha.n_cols;
   int n = fasta.size();
   int t = fasta[0].size();
   int m = t - k + 1;
   
   for (int i = 0; i < n; ++i) {
     const auto &seq = fasta[i];
     double intern_sum = 0.0;
     for (int j = 0; j < m; ++j) {
       intern_sum += probSeqGivenPosLog(seq, alpha, beta, j);
     }
     q += intern_sum;
   }
   
   return q;
 }

//'Computes de Log-Likelihood of the data and parameters.
//'@name LL
//'@param mod Can be OOPS, ZOOPS or ANR.
//'@param fasta Dataset of sequences.
//'@param alpha PWM motif model.
//'@param beta Markov backround model.
//'@param w Priori probabilities.
//'@return The Log-Likelihood of all sequences and parameters.
double LL(const std::string &mod, const std::vector<std::string> &fasta, const arma::mat &alpha, const arma::mat &beta, const double w) {
  double ll = 0.0;
  int k = alpha.n_cols;
  int n = fasta.size();
  int t = fasta[0].size();
  int m = t - k + 1;
  
  for (int i = 0; i < n; ++i) {
   const auto &seq = fasta[i];
   double intern_sum = 0.0;
   for (int j = 0; j < m; ++j) {
     if (mod == "OOPS") {
       intern_sum += probSeqGivenPos(seq, alpha, beta, j);
     }
     
     else if (mod == "ZOOPS") {
       intern_sum += w * probSeqGivenPos(seq, alpha, beta, j) + (1-w) * probSeqGivenBeta(seq, beta);
     }
     
     else if (mod == "ANR") {
       const auto &kmer = seq.substr(j, k);
       intern_sum += w * probSeqGivenAlpha(kmer, alpha) + (1-w) * probSeqGivenBeta(kmer, beta);
     }
   }
   ll += std::log(intern_sum);
  }
  
  return ll;
}

//'Create tau-order Markov Chain from the dataset.
//'@name createMarkovChain
//'@param fasta The dataset os sequences.
//'@param tau The order of markov chain.
//'@return The Markov Model of tau-order.
arma::mat createMarkovChain(const std::vector<std::string> &fasta, const int tau = 0) {
 arma::mat markov(std::pow(4, tau), 4);
 int t = fasta[0].size();
 int k = tau + 1;
 int m = t - k + 1;
 
 for (const auto &seq : fasta) {
   for (int j = 0; j < m; ++j) {
     const auto &kmer = seq.substr(j, k);
     int row = kmer2index(kmer.substr(0, k-1));
     int col = char2int(kmer[k-1]);
     markov(row, col) += 1;
   }
 }
 
 for (int i = 0; i < std::pow(4, tau); ++i) {
   markov.row(i) = markov.row(i) / arma::accu(markov.row(i));
 }
 
 return markov; 
}
 
//'Computes de probability of a sequence give the alpha model.
//'@name probSeqGivenAlpha
//'@param kmer Sequence with k size.
//'@param alpha PWM (Position Weight Matrix) model.
//'@return The probability of sequence given the current PWM model.
double probSeqGivenAlpha(const std::string &kmer, const arma::mat &alpha) {
 double prob = 1.0;
 
 int i = 0;
 for (const auto &c : kmer) {
   int idx = char2int(c);
   prob *= alpha(idx, i++);
 }
 
 return prob;
}

//'Computes de log-probability of a sequence give the alpha model.
//'@name probSeqGivenAlphaLog
//'@param kmer Sequence with k size.
//'@param alpha PWM (Position Weight Matrix) model.
//'@return The probability of sequence given the current PWM model.
double probSeqGivenAlphaLog(const std::string &kmer, const arma::mat &alpha) {
 double prob = 0.0;
 
 arma::mat alphalog = arma::log(alpha);
 
 int i = 0;
 for (const auto &c : kmer) {
   int idx = char2int(c);
   prob += alphalog(idx, i++);
 }
 
 return prob;
}

//'Computes de probability of a sequence give the beta model.
//'@name probSeqGivenBeta
//'@param seq Sequence with k size.
//'@param beta Markov Chain.
//'@return The probability of sequence given the Markov model. 
double probSeqGivenBeta(const std::string &seq, const arma::mat &beta) {
 double prob = 1.0;
 auto t = seq.size();
 
 for (int j = 0; j < t; ++j) {
    int col = 0;
    switch(seq[j]) {
      case 'A':
        col = 0;
        break;
      case 'C':
        col = 1;
        break;
      case 'G':
        col = 2;
        break;
      case 'T':
        col = 3;
        break;
    }
   prob *= beta(0, col);
 }
 
 return prob;
}

//'Computes de probability of a sequence give the beta model.
//'@name probSeqGivenBeta
//'@param seq Sequence with k size.
//'@param beta Markov Chain.
//'@return The probability of sequence given the Markov model. 
double probSeqGivenBetaLog(const std::string &seq, const arma::mat &beta) {
  double prob = 0.0;
  auto k = beta.n_cols;
  auto t = seq.size();
  auto m = t - k + 1;
  arma::mat betalog = arma::log(beta);
  
  for (int j = 0; j < m; ++j) {
    int col = 0;
    switch(seq[j]) {
    case 'A':
      col = 0;
      break;
    case 'C':
      col = 1;
      break;
    case 'G':
      col = 2;
      break;
    case 'T':
      col = 3;
      break;
    }
    prob += beta(0, col);
  }
  
  return prob;
}

//'Computes de probability of a sequence give the position.
//'@name probSeqGivenPos
//'@param seq Sequence with k size.
//'@param alpha PWM (Position Weight Matrix) model.
//'@param beta Markov Chain.
//'@param pos Position to compute probabilities.
//'@return The probability of sequence given the position. 
double probSeqGivenPos(const std::string &seq, const arma::mat &alpha, const arma::mat &beta, const int pos) {
 int k = alpha.n_cols;
 const auto &kmer = seq.substr(pos, k);
 double q = probSeqGivenBeta(seq, beta);
 return probSeqGivenAlpha(kmer, alpha) * (q/probSeqGivenBeta(kmer, beta));
}

//'Computes de log-probability of a sequence give the position.
//'@name probSeqGivenPosLog
//'@param seq Sequence with k size.
//'@param alpha PWM (Position Weight Matrix) model.
//'@param beta Markov Chain.
//'@param pos Position to compute probabilities.
//'@return The probability of sequence given the position. 

double probSeqGivenPosLog(const std::string &seq, const arma::mat &alpha, const arma::mat &beta, const int pos) {
 int k = alpha.n_cols;
 const auto &kmer = seq.substr(pos, k);
 double q = probSeqGivenBetaLog(seq, beta);
 return probSeqGivenAlphaLog(kmer, alpha) + (q - probSeqGivenBetaLog(kmer, beta));
}

//'Compute information content.
//'@name computeIC
//'@param alpha PWM model.
//'@param beta Markov 0-order model.
//'@return Information content score.
double computeIC(const arma::mat &alpha, const arma::mat &beta) {
 arma::mat alpha_log = arma::log2(alpha);
 arma::vec beta_log = arma::log2(beta.t());
 arma::mat result = alpha_log.each_col() - beta_log;
 return arma::accu(alpha % result);
}

//'Compute information content from a uniform distribution.
//'@name computeICU
//'@param alpha PWM model.
//'@return Information content score.
double computeICU(const arma::mat &alpha) {
 int k = alpha.n_cols;
 
 return 2.0*k + arma::accu(alpha % arma::log2(alpha));
}

//'Update the alpha model with all kmers in the sequence.
//'@name soft_update
//'@param alpha PWM model.
//'@param seq All kmers in the sequence.
//'@param posteriori All kmers posteriori distribution.
void update(arma::mat &alpha, const std::string &seq, const arma::rowvec &posteriori) {
 int k = alpha.n_cols;
 int t = seq.size();
 int m = t - k + 1;
 
 for (int i = 0; i < m; ++i) {
   const auto &kmer = seq.substr(i, k);
   int j = 0;
   for (const auto &c : kmer) {
     int row = 0;
     switch(c) {
     case 'A':
       row = 0;
       break;
     case 'C':
       row = 1;
       break;
     case 'G':
       row = 2;
       break;
     case 'T':
       row = 3;
       break;
     }
     alpha(row, j++) += posteriori[i];
   }
 }
}

//'Update the alpha model with all kmers in the sequence.
//'@name hard_update
//'@param alpha PWM model.
//'@param seq All kmers in the sequence.
//'@param posteriori All kmers posteriori distribution.
void hard_update(arma::mat &alpha, const std::string &seq, const arma::rowvec &posteriori) {
 int k = alpha.n_cols;
 int idx = posteriori.index_max();
 const auto &kmer = seq.substr(idx, k);
 int j = 0;
 for(const auto &c : kmer) {
   int row = char2int(c);
   alpha(row, j++) += 1;
 }
}

//'Convert kmers to PWM model.
//'@name kmers2alpha
//'@param kmers Kmers to convert.
//'@return PWM model from kmers.
arma::mat kmers2alpha(const std::vector<std::string> &kmers) {
 int k = kmers[0].size();
 int n = kmers.size();
 arma::mat alpha(4, k);
 alpha.fill(1);
 
 for (const auto &kmer : kmers) {
   for (int j = 0; j < k; ++j) {
     int c = char2int(kmer[j]);
     alpha(c, j) += 1;
   }
 }
 
 return alpha / (n + 4);
}

//'Convert consensus sequence or a simple kmer to PWM model.
//'@name consensus2alpha
//'@param consensus Consensus sequence or a simple kmer from convert to.
//'@return PWM model from consensus.
arma::mat consensus2alpha(const std::string &consensus) {
  int k = consensus.size();
  arma::mat alpha(4, k);
  alpha.fill(1);
  
  double confidence_level = 12.0;
  
  for(int i = 0; i < k; ++i) {
    char c = consensus[i];
    switch(c) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
        alpha(char2int(c), i) += confidence_level;
        
        break;
        
      case 'R': // puRine  
        // A
        alpha(0, i) += confidence_level / 2;
        // G
        alpha(2, i) += confidence_level / 2;
        
        break;
        
      case 'Y': // pYrimidine
        // C
        alpha(1, i) += confidence_level / 2;
        // T
        alpha(3, i) += confidence_level / 2;
        
        break;
        
      case 'S': //Strong interaction
        // G
        alpha(2, i) += confidence_level / 2;
        // C
        alpha(1, i) += confidence_level / 2;
        
        break;
        
      case 'W': // Weak interaction
        // A
        alpha(0, i) += confidence_level / 2;
        // T
        alpha(3, i) += confidence_level / 2;
        
        break;
        
      case 'K': // Keto
        // G
        alpha(2, i) += confidence_level / 2;
        
        // T
        alpha(3, i) += confidence_level / 2;
        
        break;
      case 'M': // aMino
        // A
        alpha(0, i) += confidence_level / 2;
        
        // C
        alpha(1, i) += confidence_level / 2;
        
        break;
        
      case 'B': // not-A, B follows A
        // C
        alpha(1, i) += confidence_level / 3;
        
        // T
        alpha(3, i) += confidence_level / 3;
        
        // G
        alpha(2, i) += confidence_level / 3;
        
        break;
        
      case 'D': // not-C, D follows C
        // G
        alpha(2, i) += confidence_level / 3;
        
        // A
        alpha(0, i) += confidence_level / 3;
        
        // T
        alpha(3, i) += confidence_level / 3;
        
        break;
        
      case 'H': // not-G, H follows G
        // A
        alpha(0, i) += confidence_level / 3;
        
        // C
        alpha(1, i) += confidence_level / 3;
        
        // T
        alpha(3, i) += confidence_level / 3;
        
        break;
        
        
      case 'V': // not-T (not-U), V follows U
        // G
        alpha(2, i) += confidence_level / 3;
        
        // C
        alpha(1, i) += confidence_level / 3;
        
        // A
        alpha(0, i) += confidence_level / 3;
        
        break;
        
      case 'N':
      case 'X':
        // A
        alpha(0, i) += confidence_level / 4;
        
        // C
        alpha(1, i) += confidence_level / 4;
        
        // G
        alpha(2, i) += confidence_level / 4;
        
        // T
        alpha(3, i) += confidence_level / 4;
        
        break;
      
      default: throw std::invalid_argument("Invalid character");
    }
  }
  
  return alpha / (confidence_level + 4);
}

//'Convert PWM model to kmers.
//'@name alpha2kmers
//'@param alpha PWM model to convert.
//'@param fasta Dataset of sequences.
//'@return The best kmers from dataset with respect to alpha.
std::vector<std::string> alpha2kmers(const arma::mat &alpha, const std::vector<std::string> &fasta) {
 int t = fasta[0].size();
 int k = alpha.n_cols;
 int m = t - k + 1;
 
 std::vector<std::string> kmers;
 arma::vec z(m);
 
 for(const auto &seq : fasta) {
   for (int j = 0; j < m; ++j) {
     const auto &kmer = seq.substr(j, k);
     z[j] = computeDKLU(alpha, kmer);
   }
   
   int idx = arma::index_max(z);
   const std::string &kmer = seq.substr(idx, k);
   kmers.push_back(kmer);
 }
 
 return kmers;
}

//'Create random initial PWM guess from fasta.
//'@name fasta2kmers
//'@param fasta Dataset of sequences.
//'@param k Size of kmers.
//'@return The best kmers from dataset with respect to alpha.
arma::mat fasta2alpha(const std::vector<std::string> &fasta, const int k) {
  int t = fasta[0].size();
  int m = t - k + 1;
  
  std::vector<std::string> kmers;
  for (const auto &seq : fasta) {
    std::uniform_int_distribution<> distrib(0, m);
    int u = distrib(gen);
    const auto &kmer = seq.substr(u, k);
    kmers.push_back(kmer);
  }
  
  arma::mat alpha = kmers2alpha(kmers);
  
  return alpha;
}


