#pragma once
#include <string>
#include <armadillo>
#include <omp.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <zlib.h>
#include <lz4.h>
#include <sstream>

void createSparseMT(const std::vector<std::string> &fasta, const int k);
void createCompactMT(const std::vector<std::string> &fasta, const int k);
