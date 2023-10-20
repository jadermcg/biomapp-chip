#pragma once

#include <armadillo>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
namespace fs = std::filesystem;

std::unique_ptr<std::vector<arma::mat>> build_models_from_sibligs(const int k);