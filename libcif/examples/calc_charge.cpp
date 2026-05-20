// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if(argc < 2) { std::cerr << "No file given!\n"; return 1; }

  auto cif = cif_data(argv[1]);

  double eps = std::atof(argv[2]);

  double charge = calc_total_charge(cif);

  if(charge <= eps) { std::cout << charge << " OK\n"; return 0; }
  else { std::cout << charge << " BAD\n"; return 1; }
}
