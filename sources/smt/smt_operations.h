#pragma once
#include <string>
#include <armadillo>
#include <stdexcept>
#include <omp.h>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include<iostream>
#include<fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

int ksearch(const std::string kmer);
tbb::concurrent_hash_map<std::string, uint64_t> khmap(const int k);
void hsib(const tbb::concurrent_hash_map <std::string, uint64_t> &hmap, const std::vector<std::string> &kmers, const int d);
std::map<std::string, std::map<std::string, int>> busca_direta(std::vector<std::string> &fasta, std::vector<std::string> &kmers, int d);
tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> kdive(const std::vector<std::string> &kmers, const int d);

