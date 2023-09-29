#include "zoops.h"
#include "utils.h"
#include "prob_utils.h"

//' Runs Expectation Maximization ZOOPS and reestimates model parameters.
//'@name zoops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability to each sequence has a motif.
//'@return Updated PWM model.
//[[Rcpp::export]]
arma::mat zoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  /**
   * Posteriori
   */
  arma::vec z(m);
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  double new_w = 0.0;
  
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    
    new_alpha.fill(1e-100);
    new_w = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      
      /**
       * E-STEP
       */
      
      std::atomic<double> marginal_prob(0.0);
      tbb::parallel_for(0, m, 1, [&](int j) { // Using Intel parallel FOR
        const auto &kmer = seq.substr(j, k);
        double local_z = w * probSeqGivenPos(seq, alpha, beta, j);
        z[j] = local_z;
        double current_value = marginal_prob.load(std::memory_order_relaxed);
        while (!marginal_prob.compare_exchange_weak(current_value, current_value + local_z, std::memory_order_relaxed));
      });
      
      // zoops
      double q = m * probSeqGivenBeta(seq, beta) * (1-w);
      
      /**
       * M-STEP
       */
      
      z = z / (marginal_prob.load() + q);
      new_w += arma::accu(z);
      update(new_alpha, seq, z.t());
    }
    double sumcol = arma::accu(new_alpha.col(0));
    alpha = new_alpha / sumcol;
    w = new_w / n;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return alpha;
}

//' Runs Expectation Maximization ZOOPS and reestimates model parameters.
//'@name logzoops
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability to each sequence has a motif.
//'@return Updated PWM model.
//[[Rcpp::export]]
arma::mat logzoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  
  /**
   * Posteriori
   */
  arma::vec z(m);
  
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  double new_w = 0.0;
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  while (true) {
    
    new_alpha.fill(1e-100);
    new_w = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      
      /**
       * E-STEP
       */
      
      std::atomic<double> marginal_prob_atomic(0.0);
      tbb::parallel_for(0, m, 1, [&](int j) {
        const auto &kmer = seq.substr(j, k);
        z[j] = std::log(w) + probSeqGivenPosLog(seq, alpha, beta, j);
        
        // Utiliza operações atômicas para atualizar marginal_prob
        double local_z_exp = std::exp(z[j]);
        double current_value = marginal_prob_atomic.load(std::memory_order_relaxed);
        while (!marginal_prob_atomic.compare_exchange_weak(current_value, current_value + local_z_exp, std::memory_order_relaxed));
      });
      
      double marginal_prob = marginal_prob_atomic.load();
      
      // zoops
      double q = std::log(m) + probSeqGivenBetaLog(seq, beta) + std::log(1-w);
      marginal_prob += std::exp(q);
      /**
       * M-STEP
       */
      z = arma::exp(z - std::log(marginal_prob));
      new_w += arma::accu(z);
      update(new_alpha, seq, z.t());
    }
    double sumcol = arma::accu(new_alpha.col(0));
    alpha = new_alpha / sumcol;
    w = new_w / n;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
  }
  
  return alpha;
}
