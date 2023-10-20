#pragma once
#include <strings.h>
#include <armadillo>
#include <vector>
#include <stdexcept>
#include<random>

double computeICU(const arma::mat &alpha);
double step_size(std::vector<double> convergence);
double computeIC(const arma::mat &alpha, const arma::mat &beta);
arma::mat kmers2alpha(const std::vector<std::string> &kmers);
arma::mat consensus2alpha(const std::string &consensus);
double probSeqGivenBeta(const std::string &seq, const arma::mat &beta);
double probSeqGivenBetaLog(const std::string &seq, const arma::mat &beta);
double probSeqGivenAlpha(const std::string &kmer, const arma::mat &alpha);
arma::mat fasta2alpha(const std::vector<std::string> &fasta, const int k);
double probSeqGivenAlphaLog(const std::string &seq, const arma::mat &alpha);
arma::mat createMarkovChain(const std::vector<std::string> &fasta, const int tau);
void update(arma::mat &alpha, const std::string &seq, const arma::rowvec &posteriori);
void hard_update(arma::mat &alpha, const std::string &seq, const arma::rowvec &posteriori);
double Q(const std::string &mod, const std::vector<std::string> &fasta, const arma::mat &alpha, const arma::mat &beta, const double w);
double LL(const std::string &mod, const std::vector<std::string> &fasta, const arma::mat &alpha, const arma::mat &beta, const double w);
std::vector<std::string> alpha2kmers(const arma::mat &alpha, const std::vector<std::string> &fasta);
double probSeqGivenPos(const std::string &seq, const arma::mat &alpha, const arma::mat &beta, const int pos);
double probSeqGivenPosLog(const std::string &seq, const arma::mat &alpha, const arma::mat &beta, const int pos);
bool hasConverged(const double cutoff, const int niter, const arma::mat &alpha, std::vector<double> &convergence, std::vector<double> &changes);