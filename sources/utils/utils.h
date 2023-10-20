#pragma once
#include <strings.h>
#include <armadillo>
#include <stdexcept>
#include <boost/filesystem.hpp>

int char2int(char c);
char int2char(int i);
double fastlog(double x);
double fastlog2(double x);
int computeSCORE(const arma::mat &alpha);
int kmer2index(const std::string &kmer);
std::string index2kmer(int index, int k);
double corr_freq(const std::string &correlation_str);
std::string corr(const std::string &a, const std::string b);
std::vector<std::string> readFasta(const std::string& filepath);
std::vector<std::string> readFasta(const std::string& filepath);
double fast_corr_freq(const std::string &a, const std::string b);
std::vector<std::string> getFilenames(const std::string& dirPath);
double computeDKLU(const arma::mat &alpha, const std::string &kmer);
int hamming_distance_parallel(const std::string &str1, const std::string &str2);
double computeDKL(const arma::mat &alpha, const arma::mat &beta, const std::string &kmer);
  
