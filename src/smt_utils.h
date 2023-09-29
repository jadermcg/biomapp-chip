#pragma once
#include <string>
#include <stdexcept>
#include <omp.h>
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>

int char2int(char c);
char int2char(int i);
uint64_t kmer2index(const std::string &kmer);
std::string index2kmer(uint64_t index, int k);
std::vector<std::string> readkmers(const std::string &path);
std::vector<std::string> readFasta(const std::string& filepath);
std::vector<std::string> getFilenames(const std::string& dirPath);
std::unordered_map<std::string, uint64_t> readhmap(const std::string &filename);