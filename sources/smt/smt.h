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
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <atomic>


void saveMT();
void processMT(const std::vector<std::string> &fasta, const int k, const int bsize);
