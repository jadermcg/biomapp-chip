#pragma once
#include <string>
#include <armadillo>
#include<iostream>
#include<fstream>
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>

void hmap(int nb, int k);
void consumer();