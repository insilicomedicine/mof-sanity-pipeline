// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"

#include <iostream>
#include <string>

int main(int argc, char **argv)
{
  if(argc < 2) { std::cerr << "No file given!\n"; return 1; }
  if(argc < 5) { std::cerr << "Usage: ./a.out filename.cif r_min r_max tolerance strategy penalty\n"; return 1; }

  auto cif = cif_data(argv[1]);

  double const param_too_close = std::atof(argv[2]);
  double const param_too_far = std::atof(argv[3]);
  double const tolerance = std::atof(argv[4]);
  std::string const strategy = argv[5];
  double const penalty = std::atof(argv[6]);


  std::cout << cif_quality(cif, param_too_close, param_too_far, tolerance, strategy, penalty) << std::endl;
}
