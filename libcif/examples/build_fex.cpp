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
  auto new_aligned = kabsch_algorithm(atoms, connection_points, ref_atoms, {0,1,2,3});
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
  if(connection_points.size() != 4) { throw std::runtime_error("FUCK YOU!!"); }

  
  double const length_short = atoms_range(new_linker.atoms[connection_points[0]], new_linker.atoms[connection_points[1]]);
  double const length_long = atoms_range(new_linker.atoms[connection_points[0]], new_linker.atoms[connection_points[2]]);

  std::string me_name = "Fe";
  double r_eff = 1.97;
  if(argc>2) { me_name = argv[2]; }
  if(argc>3) { r_eff = std::atof(argv[3]); }

  // HARDCODED ME POSITIONS
  std::vector<atom> const fract_me_positions
  {
    // LEFT a-c LAYER
    {.coord={0, 0, 0}, .name=me_name}, // real
    {.coord={0.5, 0, 0}, .name=me_name}, // real
    {.coord={1, 0, 0}, .name=me_name}, // fake
      //
    {.coord={0, 0, 0.5}, .name=me_name}, // real
    {.coord={0.5, 0, 0.5}, .name=me_name}, // real
    {.coord={1, 0, 0.5}, .name=me_name}, // fake
      //
    {.coord={0, 0, 1}, .name=me_name}, // fake
    {.coord={0.5, 0, 1}, .name=me_name}, // fake
    {.coord={1, 0, 1}, .name=me_name}, // fake
    // MIDDLE a-c LAYER
    {.coord={0, 0.5, 0}, .name=me_name}, // real
    {.coord={0.5, 0.5, 0}, .name=me_name}, // real
    {.coord={1, 0.5, 0}, .name=me_name}, // fake
      //
    {.coord={0, 0.5, 0.5}, .name=me_name}, // real
    {.coord={0.5, 0.5, 0.5}, .name=me_name}, // real
    {.coord={1, 0.5, 0.5}, .name=me_name}, // fake
      //
    {.coord={0, 0.5, 1}, .name=me_name}, // fake
    {.coord={0.5, 0.5, 1}, .name=me_name}, // fake
    {.coord={1, 0.5, 1}, .name=me_name}, // fake
    // RIGHT a-c LAYER
    {.coord={0, 1, 0}, .name=me_name}, // fake
    {.coord={0.5, 1, 0}, .name=me_name}, // fake
    {.coord={1, 1, 0}, .name=me_name}, // fake
      //
    {.coord={0, 1, 0.5}, .name=me_name}, // fake
    {.coord={0.5, 1, 0.5}, .name=me_name}, // fake
    {.coord={1, 1, 0.5}, .name=me_name}, // fake
      //
    {.coord={0, 1, 1}, .name=me_name}, // fake
    {.coord={0.5, 1, 1}, .name=me_name}, // fake
    {.coord={1, 1, 1}, .name=me_name}, // fake
      //
    // LEFT-MIDDLE a-c LAYER
    {.coord={-0.25, 0.25, -0.25}, .name=me_name}, // fake
    {.coord={0.25, 0.25, -0.25}, .name=me_name}, // fake
    {.coord={0.75, 0.25, -0.25}, .name=me_name}, // fake
    {.coord={1.25, 0.25, -0.25}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.25, 0.25}, .name=me_name}, // fake
    {.coord={0.25, 0.25, 0.25}, .name=me_name}, // real
    {.coord={0.75, 0.25, 0.25}, .name=me_name}, // real
    {.coord={1.25, 0.25, 0.25}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.25, 0.75}, .name=me_name}, // fake
    {.coord={0.25, 0.25, 0.75}, .name=me_name}, // real
    {.coord={0.75, 0.25, 0.75}, .name=me_name}, // real
    {.coord={1.25, 0.25, 0.75}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.25, 1.25}, .name=me_name}, // fake
    {.coord={0.25, 0.25, 1.25}, .name=me_name}, // fake
    {.coord={0.75, 0.25, 1.25}, .name=me_name}, // fake
    {.coord={1.25, 0.25, 1.25}, .name=me_name}, // fake
    // RIGHT-MIDDLE a-c LAYER
    {.coord={-0.25, 0.75, -0.25}, .name=me_name}, // fake
    {.coord={0.25, 0.75, -0.25}, .name=me_name}, // fake
    {.coord={0.75, 0.75, -0.25}, .name=me_name}, // fake
    {.coord={1.25, 0.75, -0.25}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.75, 0.25}, .name=me_name}, // fake
    {.coord={0.25, 0.75, 0.25}, .name=me_name}, // real
    {.coord={0.75, 0.75, 0.25}, .name=me_name}, // real
    {.coord={1.25, 0.75, 0.25}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.75, 0.75}, .name=me_name}, // fake
    {.coord={0.25, 0.75, 0.75}, .name=me_name}, // real
    {.coord={0.75, 0.75, 0.75}, .name=me_name}, // real
    {.coord={1.25, 0.75, 0.75}, .name=me_name}, // fake
      //
    {.coord={-0.25, 0.75, 1.25}, .name=me_name}, // fake
    {.coord={0.25, 0.75, 1.25}, .name=me_name}, // fake
    {.coord={0.75, 0.75, 1.25}, .name=me_name}, // fake
    {.coord={1.25, 0.75, 1.25}, .name=me_name}, // fake
  };
  
  double const a = 2.0*(length_short + r_eff) + r_eff;
  double const c = 2.0*(length_long + r_eff) + r_eff;
  double const b = 4*(std::sqrt((length_long+r_eff)*(length_long+r_eff) - 0.25*0.25*c*c)) + 2*r_eff;
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
  auto qwe0 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {1,2,4,5});
  auto qwe1 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {3,4,6,7});
  auto qwe2 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {9,10,12,13});
  auto qwe3 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {13,14,16,17});
  // VERTICAL SIDE
  auto qwe4 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {27,28,31,32});
  auto qwe5 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {32,33,36,37});
  auto qwe6 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {44,45,48,49});
  auto qwe7 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {47,48,51,52});
  // DIAGONAL SHIT
  auto qwe8 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {0,1,32,33});
  auto qwe9 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {3,4,31,32});
  auto qwe10 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {4,5,37,38});
  auto qwe11 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {7,8,36,37});
  auto qwe12 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {32,33,10,11});
  auto qwe13 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {33,34,13,14});
  auto qwe14 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {35,36,12,13});
  auto qwe15 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {36,37,15,16});
  auto qwe16 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {10,11,49,50});
  auto qwe17 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {13,14,48,49});
  auto qwe18 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {12,13,52,53});
  auto qwe19 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {15,16,51,52});
  auto qwe20 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {47,48,18,19});
  auto qwe21 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {48,49,21,22});
  auto qwe22 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {52,53,22,23});
  auto qwe23 = make_linker(new_linker.atoms, connection_points, tmp_cif.abs_atoms, {53,54,25,26});

  cif_data new_cif;
  new_cif.cell = tmp_cif.cell;
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[0]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[1]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[3]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[4]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[9]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[10]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[12]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[13]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[32]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[33]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[36]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[37]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[48]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[49]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[52]);
  new_cif.abs_atoms.push_back(tmp_cif.abs_atoms[53]);

  append_vector(qwe0, new_cif.abs_atoms);
  append_vector(qwe1, new_cif.abs_atoms);
  append_vector(qwe2, new_cif.abs_atoms);
  append_vector(qwe3, new_cif.abs_atoms);
  append_vector(qwe4, new_cif.abs_atoms);
  append_vector(qwe5, new_cif.abs_atoms);
  append_vector(qwe6, new_cif.abs_atoms);
  append_vector(qwe7, new_cif.abs_atoms);
  append_vector(qwe8, new_cif.abs_atoms);
  append_vector(qwe9, new_cif.abs_atoms);
  append_vector(qwe10, new_cif.abs_atoms);
  append_vector(qwe11, new_cif.abs_atoms);
  append_vector(qwe12, new_cif.abs_atoms);
  append_vector(qwe13, new_cif.abs_atoms);
  append_vector(qwe14, new_cif.abs_atoms);
  append_vector(qwe15, new_cif.abs_atoms);
  append_vector(qwe16, new_cif.abs_atoms);
  append_vector(qwe17, new_cif.abs_atoms);
  append_vector(qwe18, new_cif.abs_atoms);
  append_vector(qwe19, new_cif.abs_atoms);
  append_vector(qwe20, new_cif.abs_atoms);
  append_vector(qwe21, new_cif.abs_atoms);
  append_vector(qwe22, new_cif.abs_atoms);
  append_vector(qwe23, new_cif.abs_atoms);

  new_cif.fract_atoms.clear();
  new_cif.calc_fract_coords();

  std::ofstream out("xef.cif");
  print_cif(new_cif, &out);
}
