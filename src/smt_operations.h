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

tbb::concurrent_hash_map<std::string, uint64_t> khmap(const int k);
tbb::concurrent_hash_map<std::string, uint64_t> hmap();
tbb::concurrent_hash_map<std::string, tbb::concurrent_hash_map<std::string, uint64_t>> kdive(const std::vector<std::string> &kmers, const int d);
void hsib(const tbb::concurrent_hash_map <std::string, uint64_t> &hmap, const std::vector<std::string> &kmers, const int d);

