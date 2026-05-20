// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "graph_ops.hpp"
#include "transform_ops.hpp"
#include "utils.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<atom> make_linker(std::vector<atom> const& atoms, std::vector<std::size_t> const& connection_points,
                              std::vector<atom> const& abs_atoms,
                              std::vector<std::size_t> const& ref_me_indices,
                              std::vector<vec3d> const& shifts = {})
{
  if(shifts.size() != ref_me_indices.size() && shifts.size() != 0) { throw std::runtime_error("wrong shifts!"); }

  std::vector<atom> ref_atoms;
  for(std::size_t index : ref_me_indices) { ref_atoms.push_back(abs_atoms[index]); }
  for(std::size_t i=0; i<shifts.size(); ++i) { ref_atoms[i].coord += shifts[i]; }
  auto new_aligned = kabsch_algorithm(atoms, connection_points, ref_atoms, {0,1,2});
  return new_aligned;
}

int main(int argc, char** argv)
{
  if(argc < 2) { return 1; }

  std::ifstream ifs(argv[1]);
  auto [atoms, comment] = read_xyz(ifs);
  auto new_linker = linker {.atoms=atoms};
  new_linker.graph = calc_connect_graph_v2(new_linker.atoms, new_linker.atoms.size(), 1.3, "cordero");

  std::vector<std::size_t> connection_points; // IN ORDER!!!
  std::stringstream ss(comment);
  int tmp;
  while(ss>>tmp) { connection_points.push_back(tmp); }
  if(connection_points.size() != 3) { throw std::runtime_error("FUCK YOU!!"); }

  
  double const length_short = atoms_range(new_linker.atoms[connection_points[1]], new_linker.atoms[connection_points[2]]);
  double const length_long = atoms_range(new_linker.atoms[connection_points[0]], new_linker.atoms[connection_points[2]]);

  std::string me_name = "Ni";
  double r_eff = 1.97;
  if(argc>2) { me_name = argv[2]; }
  if(argc>3) { r_eff = std::atof(argv[3]); }

  // HARDCODED ME POSITIONS
  std::vector<atom> const fract_me_positions
  {
    // Left
    {.coord={0,0,0}, .name=me_name}, // real
    {.coord={1,0,0}, .name=me_name}, // fake
    {.coord={0,1,0}, .name=me_name}, // fake
    {.coord={1,1,0}, .name=me_name}, // fake
    // Right
    {.coord={0,0,1}, .name=me_name}, // fake
    {.coord={1,0,1}, .name=me_name}, // fake
    {.coord={0,1,1}, .name=me_name}, // fake
    {.coord={1,1,1}, .name=me_name}, // fake
    // Middle
    {.coord={-0.5,0.5,0.5}, .name=me_name}, // fake
    {.coord={0.5,0.5,0.5}, .name=me_name}, // real
    {.coord={1.5,0.5,0.5}, .name=me_name}, // fake
  };
  
  double const a = length_short + 1.5*r_eff;
  double const c = 2.0*(length_long) + 1.5*r_eff;
  double const b = 2*(std::sqrt((length_long+r_eff)*(length_long+r_eff) - 0.5*0.5*c*c)) + 2*r_eff;
  std::cerr << "a = " << a << std::endl;
  std::cerr << "b = " << b << std::endl;
  std::cerr << "c = " << c << std::endl;


  double const alpha = 90;
  double const beta = 90;
  double const gamma = 90;

  cif_data tmp_cif;
  tmp_cif.cell = make_unit_cell(a,b,c,alpha,beta,gamma);
  append_vector(fract_me_positions, tmp_cif.fract_atoms);
  tmp_cif.calc_abs_coords();

  // VERTICAL MAIN
  auto qwe0 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {9,0,1});
  auto qwe1 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {2,9,8});
  auto qwe2 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {4,8,9});
  auto qwe3 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {9,7,6});

  cif_data new_cif;
  new_cif.cell = tmp_cif.cell;
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[0]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[9]);

  append_vector(qwe0, new_cif.abs_atoms);
  append_vector(qwe1, new_cif.abs_atoms);
  append_vector(qwe2, new_cif.abs_atoms);
  append_vector(qwe3, new_cif.abs_atoms);

  new_cif.fract_atoms.clear();
  new_cif.calc_fract_coords();

  std::ofstream out("tub.cif");
  print_cif(new_cif, &out);
}
