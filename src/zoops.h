#pragma once
#include <RcppArmadillo.h>
#include <strings.h>
#include <vector>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <atomic>
//[[Rcpp::depends(RcppArmadillo)]]

arma::mat logzoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
arma::mat zoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
