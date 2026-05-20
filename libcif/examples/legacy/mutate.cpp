// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"
#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
  if(argc < 2) { std::cout << "No file given!\n"; return 1; }
  // if(argc < 5) { std::cerr << "Usage: ./a.out filename.cif r_min r_max tolerance\n"; return 1; }

  std::ifstream inp(argv[1]);
  auto cif = read_mutant(inp);


  std::vector<int> new_mutations;
  for(int i=2; i<argc; ++i)
  {
    new_mutations.push_back(std::stoi(argv[i]));
  }

  if(new_mutations.size() != cif.cmsm.genome.size()) { std::cout << "Wrong chromosome size!!\n"; return 1; }

  mutate(cif, new_mutations);
  // cif.calc_fract_coords();
  print_cif(cif);
}
