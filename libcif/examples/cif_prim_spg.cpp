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
  double tolerance = 0.25;
  if(argc > 2) { tolerance = std::atof(argv[2]); }

  cell_prim_spg(cif, tolerance);
  print_cif(cif);
}
