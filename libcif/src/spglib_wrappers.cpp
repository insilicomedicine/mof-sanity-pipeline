// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "spglib_wrappers.hpp"

#include "cif.hpp"
#include "spglib.h"
#include <iostream>
#include <unordered_map>


namespace spglib_wrappers
{

std::array<std::array<double,3>,3> make_axes_matrix_sasha(unit_cell const& cell)
{
  double rad_alpha = cell.alpha * M_PI / 180;
  double rad_beta = cell.beta * M_PI / 180;
  double rad_gamma = cell.gamma * M_PI / 180;

  double v = std::sqrt(1.0 - std::cos(rad_alpha) * std::cos(rad_alpha) -
                       std::cos(rad_beta) * std::cos(rad_beta) -
                       std::cos(rad_gamma) * std::cos(rad_gamma) +
                       2.0 * std::cos(rad_alpha) * std::cos(rad_beta) * std::cos(rad_gamma));

  std::array<std::array<double,3>,3> lattice;
  lattice[0][0] = cell.a;
  lattice[1][0] = cell.b * std::cos(rad_gamma);
  lattice[2][0] = cell.c * std::cos(rad_beta),
  
  lattice[0][1] = 0;
  lattice[1][1] = cell.b * std::sin(rad_gamma);
  lattice[2][1] = cell.c * (std::cos(rad_alpha) - std::cos(rad_beta) * std::cos(rad_gamma)) / std::sin(rad_gamma);
  
  lattice[0][2] = 0;
  lattice[1][2] = 0;
  lattice[2][2] = cell.c * v / std::sin(rad_gamma);

  return lattice;
}

void primitivize_cell(cif_data& cif, double tolerance)
{
  auto lattice_matrix = make_axes_matrix_sasha(cif.cell);
  double c_lat_mat[3][3];
  for(int k=0; k<3; ++k)
  {
    for(int l=0; l<3; ++l)
    {
      c_lat_mat[k][l] = lattice_matrix[l][k];
    }
  }
  
  std::vector<double[3]> positions(cif.num_atoms);
  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    for(int k=0; k<3; ++k)
    {
      positions[i][k] = cif.fract_atoms[i].coord[k];
    }
  }

  int curr_type = 0;
  std::unordered_map<std::string,int> types_map;
  for(auto const& atom : cif.fract_atoms)
  {
    if(!types_map.contains(atom.name)) { types_map[atom.name] = curr_type++; }
  }

  std::vector<int> types(cif.fract_atoms.size());
  for(std::size_t i=0; i<cif.fract_atoms.size(); ++i)
  {
    types[i] = types_map[cif.fract_atoms[i].name];
  }

  auto new_num_atoms = spg_standardize_cell(c_lat_mat, positions.data(), types.data(), cif.num_atoms, true, false, tolerance);
  spg_niggli_reduce(c_lat_mat, tolerance);

  std::cerr << "old " << cif.num_atoms << std::endl;
  std::cerr << "new_num_atoms " << new_num_atoms << std::endl;
  
  for(int k=0; k<3; ++k)
  {
    for(int l=0; l<3; ++l)
    {
      lattice_matrix[k][l] = c_lat_mat[l][k];
      std::cerr << lattice_matrix[k][l] << ' ';
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  crystal_axes new_axes {.a=lattice_matrix[0], .b=lattice_matrix[1], .c=lattice_matrix[2]};
  auto new_cell = make_unit_cell(new_axes);
  cif.cell = new_cell;

  std::vector<atom> new_atoms(new_num_atoms);
  for(int i=0; i<new_num_atoms; ++i)
  {
    for(int k=0; k<3; ++k)
    {
      new_atoms[i].coord[k] = positions[i][k];
    }

    for(auto const& [key, value] : types_map)
    {
      if(value == types[i]) { new_atoms[i].name = key; }
    }
  }

  cif.num_atoms = new_num_atoms;
  cif.fract_atoms = new_atoms;
}

}
