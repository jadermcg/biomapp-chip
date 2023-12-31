#pragma once
#include <armadillo>
#include <string.h>
#include <vector>
#include <omp.h>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <atomic>

arma::mat oops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w);
arma::mat logoops(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, double cutoff, int niter, double w);
