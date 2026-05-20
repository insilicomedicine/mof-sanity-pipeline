// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if(argc < 2) { std::cerr << "No file given!\n"; return 1; }

  double const tolerance = std::atof(argv[1]);
  std::string strategy {argv[2]};

  for(int n=3; n<argc; ++n)
  {
    auto cif = cif_data(argv[n]);


    std::cout << argv[n] << std::endl;
    std::cout << calc_cif_hash(cif, false, true, tolerance, strategy) << '\n';
    std::cout << calc_cif_hash(cif, false, false, tolerance, strategy) <<'\n';
    std::cout << calc_cif_hash(cif, true, true, tolerance, strategy) << '\n';
    std::cout << calc_cif_hash(cif, true, false, tolerance, strategy) << '\n';
    std::cout << '\n';
  }
}
