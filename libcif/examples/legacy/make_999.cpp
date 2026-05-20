// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "build_ops.hpp"
#include "cif.hpp"
#include <cmath>

int main(int argc, char** argv)
// int main()
{
  if(argc < 3) { return 1; }
  auto node = read_xyz_to_node(argv[1]);
  auto edge = read_xyz_to_edge(argv[2]);

  constexpr double const third = 1.0/3.0;

  auto topology = building_blocks::topology
  {
    .nodes =
    {
      {.type = node, .frac_coord = {0,0,0}, .orientation = {0,0,0}},
      {.type = node, .frac_coord = {third,third,0}, .orientation = {M_PI,0,0}},

      {.type = node, .frac_coord = {third,third-1,0}, .orientation = {M_PI,0,0}},
      {.type = node, .frac_coord = {third-1,third,0}, .orientation = {M_PI,0,0}},

      {.type = node, .frac_coord = {1,0,0}, .orientation = {0,0,0}},
      {.type = node, .frac_coord = {0,1,0}, .orientation = {0,0,0}},
    },
    .edges = 
    {
      {.type=edge, .connections{{0,0}, {1,0}}, .normal_to={0,0,1}},

      {.type=edge, .connections{{0,1}, {2,1}}, .normal_to={0,0,1}},
      {.type=edge, .connections{{0,2}, {3,2}}, .normal_to={0,0,1}},

      {.type=edge, .connections{{1,1}, {5,1}}, .normal_to={0,0,1}},
      {.type=edge, .connections{{1,2}, {4,2}}, .normal_to={0,0,1}},

    }
  };

  auto cif = build_999(topology);
  print_cif(cif);
}
