// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "cif_ops.hpp"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if(argc < 2) { std::cerr << "No file given!\n"; return 1; }
  if(argc < 5) { std::cerr << "Usage: ./a.out filename.cif r_min r_max tolerance\n"; return 1; }

  auto cif = cif_data(argv[1]);

  double const param_too_close = std::atof(argv[2]);
  double const param_too_far = std::atof(argv[3]);
  double const tolerance = std::atof(argv[4]);
  
  std::string strategy {"cordero"};
  if(argc > 5) { strategy = argv[5]; }


  std::string formula = calc_cif_formula(cif);
  double V = calc_volume(cif);
  double rho = calc_density(cif);
  bool good = cif_is_good(cif, param_too_close, param_too_far, tolerance, strategy);

  std::cout << (good ? "OK " : "BAD ");
  std::cout << formula << ' ' << V << ' ' << rho << ' ';
  std::cout << cif.cell.a << ' ' << cif.cell.b << ' ' << cif.cell.c << ' ';
  std::cout << cif.cell.alpha << ' ' << cif.cell.beta << ' ' << cif.cell.gamma << std::endl;
  
  return good ? 0 : 1;
}
