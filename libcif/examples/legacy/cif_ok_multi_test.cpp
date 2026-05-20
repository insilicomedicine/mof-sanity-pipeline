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


  double const param_too_close = std::atof(argv[1]);
  double const param_too_far = std::atof(argv[2]);
  double const tolerance = std::atof(argv[3]);

  std::string strategy {"cordero"};
  // if(argc > 5) { strategy = argv[5]; }

  for(int n=4; n<argc; ++n)
  {
    auto cif = cif_data(argv[n]);
    bool good = cif_is_good(cif, param_too_close, param_too_far, tolerance, strategy);
    if(good)
    {
      std::cout << argv[n] << " OK\n";
    }
    else
    {
      std::cout << argv[n] << " BAD\n";
    }
  }
}
