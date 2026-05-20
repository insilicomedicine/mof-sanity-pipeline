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

  auto cif = cif_data(argv[1]);

  double const tolerance = std::atof(argv[2]);
  
  std::string strategy {"cordero"};
  if(argc > 3) { strategy = argv[3]; }


  std::cout << calc_cif_hash(cif, false, true, tolerance, strategy) << std::endl;
  std::cout << calc_cif_hash(cif, false, false, tolerance, strategy) << std::endl;
  std::cout << calc_cif_hash(cif, true, true, tolerance, strategy) << std::endl;
  std::cout << calc_cif_hash(cif, true, false, tolerance, strategy) << std::endl;
}
