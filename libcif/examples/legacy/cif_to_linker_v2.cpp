// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"

#include <iostream>
#include <string>

int main(int argc, char **argv)
{
  if(argc < 2) { std::cout << "No file given!\n"; return 1; }
  if(argc < 5) { std::cerr << "Usage: ./a.out filename.cif r_min r_max tolerance\n"; return 1; }

  auto cif = cif_data(argv[1]);

  double const param_too_close = std::atof(argv[2]);
  double const param_too_far = std::atof(argv[3]);
  double const tolerance = std::atof(argv[4]);
  
  std::string strategy {"cordero"};
  if(argc > 5) { strategy = argv[5]; }

  auto uniq_linkers = separate_linkers_smiles_v2(cif, 0, param_too_close, param_too_far, tolerance, strategy);

  for(auto const& uniq_linker : uniq_linkers)
  {
      std::cout << uniq_linker << std::endl;
  }
}
