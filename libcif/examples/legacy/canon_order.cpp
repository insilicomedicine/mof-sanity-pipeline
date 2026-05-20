// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"
#include "utils.hpp"
#include "rdkit_wrappers.hpp"

#include <iostream>
#include <string>

std::vector<atom> sort_by_canon_order(std::vector<atom> const& atoms, std::vector<std::size_t> canon_order)
{
  std::map<std::size_t, atom> new_order;
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    new_order[canon_order[i]] = atoms[i];
  }

  std::vector<atom> new_atoms;
  for(auto const& [key,value] : new_order)
  {
    new_atoms.push_back(value);
  }
  return new_atoms;
}

int main(int argc, char **argv)
{
  if(argc < 2) { std::cout << "No file given!\n"; return 1; }
  if(argc < 5) { std::cerr << "Usage: ./a.out filename.cif r_min r_max tolerance\n"; return 1; }

  auto cif = cif_data(argv[1]);

  double const param_too_close = std::atof(argv[2]);
  double const param_too_far = std::atof(argv[3]);
  double const tolerance = std::atof(argv[4]);
  
  // int const charge = std::atoi(argv[5]);
  std::string strategy {"cordero"};
  // if(argc > 5) { strategy = argv[5]; }

  auto all_linkers = separate_linkers(cif, param_too_close, param_too_far, tolerance, strategy);

  for(auto& linker : all_linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    auto canon_order = rdkit_wrappers::calc_canon_order_rdkit(linker.atoms);
    print_xyz(sort_by_canon_order(linker.atoms, canon_order));
  }

  // for(auto& linker: all_linkers)
  // {
  //   linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
  //   print_xyz(linker.atoms);
  //   auto canon_order = calc_canon_order_rdkit(linker.atoms);
  //   for(auto i : canon_order)
  //   {
  //       std::cout << i << ' ';
  //   }
  //   std::cout << std::endl;
  // }
}
