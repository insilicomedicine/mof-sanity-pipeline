// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "build_ops.hpp"
#include "cif.hpp"

int main(int argc, char** argv)
{
  if(argc < 2) { return 1; }
  auto edge = read_xyz_to_edge(argv[1]);

  auto topology = building_blocks::topology
  {
    .nodes =
    {
      {.type = 0, .frac_coord = {0,0,0}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {1,0,0}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {-1,0,0}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {0,1,0}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {0,-1,0}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {0,0,1}, .orientation = {0,0,0}},
      {.type = 0, .frac_coord = {0,0,-1}, .orientation = {0,0,0}}
    },
    .edges = 
    {
      {.type=edge, .connections{{0,0}, {1,1}}, .normal_to={0,1,0}},
      {.type=edge, .connections{{0,1}, {2,0}}, .normal_to={0,1,0}},
      {.type=edge, .connections{{0,2}, {3,3}}, .normal_to={0,0,1}},
      {.type=edge, .connections{{0,3}, {4,2}}, .normal_to={0,0,1}},
      {.type=2, .connections{{0,4}, {5,5}}, .normal_to={0,1,0}},
      {.type=2, .connections{{0,5}, {6,4}}, .normal_to={0,1,0}},
    }
    // .edges = 
    // {
    //   {.type=0, .connections{{0,0}, {1,1}}, .orientation={0,1,0}},
    //   {.type=0, .connections{{0,1}, {2,0}}, .orientation={0,1,0}},
    //   {.type=0, .connections{{0,2}, {3,3}}, .orientation={0,0,1}},
    //   {.type=0, .connections{{0,3}, {4,2}}, .orientation={0,0,1}},
    //   {.type=2, .connections{{0,4}, {5,5}}, .orientation={0,1,0}},
    //   {.type=2, .connections{{0,5}, {6,4}}, .orientation={0,1,0}},
    // }
  };

  auto cif = build_hum(topology);
  print_cif(cif);
}
