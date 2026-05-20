// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "graph_ops.hpp"
#include "transform_ops.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>
#include <string>

std::vector<atom> old_lower_node
{
{ .coord={ 1.493980000, 3.160920000, 1.821050000 }, .name="H"},
{ .coord={ -3.174170000, 1.086690000, 4.420850000 }, .name="H"},
{ .coord={ 3.342080000, 3.174320000, 3.812600000 }, .name="H"},
{ .coord={ 8.004790000, 1.092560000, 1.196380000 }, .name="H"},
{ .coord={ 3.046210000, 1.819780000, 1.196920000 }, .name="H"},
{ .coord={ 1.790290000, 1.803110000, 4.423850000 }, .name="H"},
{ .coord={ 1.198360000, 2.266640000, 1.586240000 }, .name="O"},
{ .coord={ 5.767620000, 1.118420000, 4.429090000 }, .name="O"},
{ .coord={ -4.077990000, 3.816150000, 4.793640000 }, .name="O"},
{ .coord={ 5.850770000, 3.163280000, 2.621730000 }, .name="O"},
{ .coord={ 3.639280000, 2.278130000, 4.039830000 }, .name="O"},
{ .coord={ -1.327810000, 0.627833000, 4.036390000 }, .name="O"},
{ .coord={ 6.153360000, 0.631021000, 1.583290000 }, .name="O"},
{ .coord={ -0.936795000, 1.123720000, 1.189560000 }, .name="O"},
{ .coord={ 4.027430000, 1.791030000, 1.189940000 }, .name="O"},
{ .coord={ 0.809485000, 1.770460000, 4.430130000 }, .name="O"},
{ .coord={ 8.928500000, 3.827340000, 0.846995000 }, .name="O"},
{ .coord={ -1.003950000, 3.156080000, 3.013010000 }, .name="O"},
{ .coord={ -0.005723390, 2.325270000, 0.085763400 }, .name="Al"},
{ .coord={ 4.896960000, 1.456970000, 2.810330000 }, .name="Al"},
{ .coord={ 4.949470000, 0.585116000, 0.082653700 }, .name="Al"},
{ .coord={ -0.064255000, 1.444690000, 2.810550000 }, .name="Al"},
{ .coord={ -0.994209000, 3.827340000, 0.846995000 }, .name="O"},
{ .coord={ -4.071940000, 3.163280000, 2.621730000 }, .name="O"},
{ .coord={ 8.918760000, 3.156080000, 3.013010000 }, .name="O"},
{ .coord={ 5.844720000, 3.816150000, 4.793640000 }, .name="O"}
};

std::vector<atom> old_linker0
{
{ .coord={ -0.331052000, 6.432690000, 1.680080000 }, .name="H"},
{ .coord={ 0.725972000, 8.010540000, 1.852530000 }, .name="H"},
{ .coord={ 0.836974000, 9.508470000, 2.777820000 }, .name="H"},
{ .coord={ 0.485859000, 9.567140000, 1.053490000 }, .name="H"},
{ .coord={ -0.339930000, 11.449200000, 1.700430000 }, .name="H"},
{ .coord={ -4.077990000, 3.816150000, 4.793640000 }, .name="O"},
{ .coord={ -1.003950000, 3.156080000, 3.013010000 }, .name="O"},
{ .coord={ -4.132970000, 14.049000000, 4.747160000 }, .name="O"},
{ .coord={ -1.045090000, 14.725300000, 2.983120000 }, .name="O"},
{ .coord={ -1.861750000, 7.727330000, 2.546390000 }, .name="C"},
{ .coord={ -1.288450000, 6.465440000, 2.197440000 }, .name="C"},
{ .coord={ -1.916070000, 5.273560000, 2.470060000 }, .name="C"},
{ .coord={ -1.888140000, 10.160900000, 2.535260000 }, .name="C"},
{ .coord={ -1.184280000, 8.950270000, 2.328650000 }, .name="C"},
{ .coord={ 0.280946000, 9.001740000, 1.972990000 }, .name="C"},
{ .coord={ -1.313870000, 11.420300000, 2.187730000 }, .name="C"},
{ .coord={ -1.957710000, 12.609800000, 2.437760000 }, .name="C"},
{ .coord={ -1.268510000, 3.983570000, 2.081400000 }, .name="C"},
{ .coord={ -3.801420000, 3.980390000, 3.561040000 }, .name="C"},
{ .coord={ -0.994209000, 3.827340000, 0.846995000 }, .name="O"},
{ .coord={ -1.041640000, 14.060800000, 0.814807000 }, .name="O"},
{ .coord={ -3.199280000, 12.605200000, 3.132600000 }, .name="C"},
{ .coord={ -3.822680000, 11.409400000, 3.405490000 }, .name="C"},
{ .coord={ -4.789010000, 11.430200000, 3.908510000 }, .name="H"},
{ .coord={ -3.233170000, 10.154100000, 3.067680000 }, .name="C"},
{ .coord={ -3.909920000, 8.932610000, 3.293930000 }, .name="C"},
{ .coord={ -3.208130000, 7.722600000, 3.091360000 }, .name="C"},
{ .coord={ -3.776900000, 6.462900000, 3.453340000 }, .name="C"},
{ .coord={ -4.736770000, 6.437250000, 3.967430000 }, .name="H"},
{ .coord={ -3.148010000, 5.270960000, 3.183120000 }, .name="C"},
{ .coord={ -1.313120000, 13.901900000, 2.049160000 }, .name="C"},
{ .coord={ -4.071940000, 3.163280000, 2.621730000 }, .name="O"},
{ .coord={ -5.366710000, 8.929310000, 3.685530000 }, .name="C"},
{ .coord={ -5.866730000, 8.012420000, 3.353610000 }, .name="H"},
{ .coord={ -5.896990000, 9.768360000, 3.217220000 }, .name="H"},
{ .coord={ -5.528810000, 9.010300000, 4.770180000 }, .name="H"},
{ .coord={ -3.854270000, 13.893700000, 3.513680000 }, .name="C"},
{ .coord={ -4.117960000, 14.719700000, 2.580810000 }, .name="O"},
};

std::vector<atom> old_middle_node
{
{ .coord={ 3.304910000, 14.717500000, 3.775880000 }, .name="H"},
{ .coord={ 1.461700000, 14.704700000, 1.780790000 }, .name="H"},
{ .coord={ 6.429690000, 18.144000000, 1.787230000 }, .name="H"},
{ .coord={ -1.647290000, 18.148100000, 3.782320000 }, .name="H"},
{ .coord={ 7.976980000, 16.793200000, 1.173350000 }, .name="H"},
{ .coord={ -3.199830000, 16.784100000, 4.400120000 }, .name="H"},
{ .coord={ 1.757010000, 16.061900000, 4.398980000 }, .name="H"},
{ .coord={ 3.011280000, 16.076800000, 1.172750000 }, .name="H"},
{ .coord={ 3.603410000, 15.610900000, 4.011130000 }, .name="O"},
{ .coord={ 1.163680000, 15.600300000, 1.554760000 }, .name="O"},
{ .coord={ 6.132020000, 17.249100000, 1.558040000 }, .name="O"},
{ .coord={ -1.350590000, 17.253400000, 4.014350000 }, .name="O"},
{ .coord={ -0.964684000, 16.760900000, 1.166150000 }, .name="O"},
{ .coord={ 5.742130000, 16.749600000, 4.406340000 }, .name="O"},
{ .coord={ 0.775825000, 16.090100000, 4.406560000 }, .name="O"},
{ .coord={ 3.992130000, 16.111000000, 1.166540000 }, .name="O"},
{ .coord={ 8.881070000, 14.060800000, 0.814807000 }, .name="O"},
{ .coord={ -4.132970000, 14.049000000, 4.747160000 }, .name="O"},
{ .coord={ 0.852630000, 18.786600000, 4.756490000 }, .name="O"},
{ .coord={ 3.944490000, 18.816400000, 0.823372000 }, .name="O"},
{ .coord={ -1.045090000, 14.725300000, 2.983120000 }, .name="O"},
{ .coord={ 5.804750000, 14.719700000, 2.580810000 }, .name="O"},
{ .coord={ 0.863631000, 18.123900000, 2.587090000 }, .name="O"},
{ .coord={ 3.936190000, 18.143900000, 2.989440000 }, .name="O"},
{ .coord={ -0.042110600, 15.560400000, 0.055811700 }, .name="Al"},
{ .coord={ 4.926590000, 17.303500000, 0.058921400 }, .name="Al"},
{ .coord={ -0.093056100, 16.425000000, 2.786000000 }, .name="Al"},
{ .coord={ 4.868630000, 16.431300000, 2.785400000 }, .name="Al"},
{ .coord={ -1.041640000, 14.060800000, 0.814807000 }, .name="O"},
{ .coord={ -4.117960000, 14.719700000, 2.580810000 }, .name="O"},
{ .coord={ 5.789740000, 14.049000000, 4.747160000 }, .name="O"},
{ .coord={ 8.877620000, 14.725300000, 2.983120000 }, .name="O"},
};

std::vector<atom> old_linker1
{
{ .coord={ 0.115796000, 26.404100000, 3.979650000 }, .name="H"},
{ .coord={ 4.513300000, 26.436700000, 1.676530000 }, .name="H"},
{ .coord={ -0.912243000, 24.810900000, 3.869770000 }, .name="H"},
{ .coord={ -0.632242000, 23.252800000, 4.649940000 }, .name="H"},
{ .coord={ -1.038220000, 23.317400000, 2.938600000 }, .name="H"},
{ .coord={ 5.658130000, 24.884000000, 2.293790000 }, .name="H"},
{ .coord={ 5.362840000, 23.824700000, 0.912736000 }, .name="H"},
{ .coord={ 5.709660000, 23.135100000, 2.501710000 }, .name="H"},
{ .coord={ 4.606190000, 21.452200000, 1.704080000 }, .name="H"},
{ .coord={ 0.160456000, 21.379900000, 3.919470000 }, .name="H"},
{ .coord={ 3.824510000, 29.041600000, 0.836466000 }, .name="O"},
{ .coord={ 0.740590000, 29.021500000, 4.784470000 }, .name="O"},
{ .coord={ 3.820020000, 29.703500000, 3.005700000 }, .name="O"},
{ .coord={ 0.752981000, 29.682500000, 2.614800000 }, .name="O"},
{ .coord={ 0.852630000, 18.786600000, 4.756490000 }, .name="O"},
{ .coord={ 3.944490000, 18.816400000, 0.823372000 }, .name="O"},
{ .coord={ 0.863631000, 18.123900000, 2.587090000 }, .name="O"},
{ .coord={ 3.936190000, 18.143900000, 2.989440000 }, .name="O"},
{ .coord={ 1.072080000, 26.377900000, 3.459610000 }, .name="C"},
{ .coord={ 3.001870000, 25.137600000, 2.563300000 }, .name="C"},
{ .coord={ 3.556010000, 26.401500000, 2.194930000 }, .name="C"},
{ .coord={ 2.912620000, 27.586800000, 2.457510000 }, .name="C"},
{ .coord={ 1.682690000, 27.575500000, 3.173900000 }, .name="C"},
{ .coord={ 1.657700000, 25.120800000, 3.113620000 }, .name="C"},
{ .coord={ 0.993001000, 23.892200000, 3.334190000 }, .name="C"},
{ .coord={ -0.461437000, 23.825100000, 3.727820000 }, .name="C"},
{ .coord={ 1.699890000, 22.688200000, 3.099760000 }, .name="C"},
{ .coord={ 3.042760000, 22.707000000, 2.563520000 }, .name="C"},
{ .coord={ 3.714850000, 23.934400000, 2.357950000 }, .name="C"},
{ .coord={ 5.177520000, 23.947600000, 1.989530000 }, .name="C"},
{ .coord={ 3.637570000, 21.459500000, 2.203330000 }, .name="C"},
{ .coord={ 3.015410000, 20.257500000, 2.450910000 }, .name="C"},
{ .coord={ 1.774750000, 20.240000000, 3.147280000 }, .name="C"},
{ .coord={ 1.130550000, 21.422300000, 3.426010000 }, .name="C"},
{ .coord={ 3.670340000, 18.971400000, 2.057830000 }, .name="C"},
{ .coord={ 1.129740000, 18.945500000, 3.523010000 }, .name="C"},
{ .coord={ 3.552410000, 28.881000000, 2.070650000 }, .name="C"},
{ .coord={ 1.021020000, 28.861600000, 3.551820000 }, .name="C"}
};

std::vector<atom> old_upper_node
{
{ .coord={ -1.747210000, 29.683400000, 3.806920000 }, .name="H"},
{ .coord={ 6.329350000, 29.688600000, 1.815980000 }, .name="H"},
{ .coord={ 3.824510000, 29.041600000, 0.836466000 }, .name="O"},
{ .coord={ 0.740590000, 29.021500000, 4.784470000 }, .name="O"},
{ .coord={ 3.820020000, 29.703500000, 3.005700000 }, .name="O"},
{ .coord={ 0.752981000, 29.682500000, 2.614800000 }, .name="O"},
};

void order_coo_groups(std::vector<std::size_t>& target, std::vector<atom> const& atoms)
{
  std::cerr << "Old order: ";
  for(std::size_t item : target) { std::cerr << item << ' '; }
  std::cerr << '\n';
  std::sort(std::next(std::begin(target)), std::end(target),
            [&atoms, &target](std::size_t lhs, std::size_t rhs)
                  { return atoms_range(atoms[lhs], atoms[target[0]]) <
                          atoms_range(atoms[rhs], atoms[target[0]]);});
  std::cerr << "New order: ";
  for(std::size_t item : target) { std::cerr << item << ' '; }
  std::cerr << '\n';
}

std::vector<std::size_t> find_coo_groups(std::vector<atom> const& atoms, molecular_graph const& graph)
{
  std::vector<std::size_t> coo_groups;
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    if(is_COO_group(i, atoms, graph)) { coo_groups.push_back(i); }
  }
  return coo_groups;
}

std::vector<atom> delete_o_from_coo(std::vector<atom> const& atoms, molecular_graph const& graph, std::vector<std::size_t> const& targets)
{
  std::unordered_set<std::size_t> os;
  for(std::size_t target : targets)
  {
    for(auto const& bond : graph.adjacency_list[target])
    {
      if(atoms[bond.neighbour].name == "O") { os.insert(bond.neighbour); }
    }
  }

  std::cerr << "Number of deleted oxygens " << os.size() << std::endl;

  std::vector<atom> new_atoms;
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    if(os.contains(i)) { continue; }
    new_atoms.push_back(atoms[i]);
  }

  std::cerr << atoms.size() << ' ' << new_atoms.size() << std::endl;
  return new_atoms;
}

double min_y_coordinate(std::vector<atom> const& atoms)
{
  auto min_atom = *std::min_element(std::begin(atoms), std::end(atoms),
                                    [](atom const& a, atom const& b){ return a.coord[1] < b.coord[1]; });
  return min_atom.coord[1];
}

int main(int argc, char** argv)
{
  if(argc < 2) { return 1; }

  // READ NEW LINKER
  std::ifstream ifs(argv[1]); // new linker
  auto [new_atoms, comment] = read_xyz(ifs);
  auto new_atoms_graph = calc_connect_graph_v2(new_atoms, new_atoms.size(), 1.3, "cordero");
  linker new_linker { .atoms = new_atoms, .graph=new_atoms_graph };
  auto new_coo_groups = find_coo_groups(new_linker.atoms, new_linker.graph);
  order_coo_groups(new_coo_groups, new_linker.atoms);

  // As per new numeration, we name new len between first and third atom
  double new_length = atoms_range(new_atoms[new_coo_groups[0]], new_atoms[new_coo_groups[2]]);

  std::vector<linker> old_linkers { linker{.atoms=old_linker0},  linker{.atoms=old_linker1}};
  for(auto& old_linker : old_linkers)
  {
    old_linker.graph = calc_connect_graph_v2(old_linker.atoms, old_linker.atoms.size(), 1.3, "cordero");
  }

  // ACTUAL ZOMBIE TRANSFORMATION
  std::vector<linker> new_linkers;
  for(auto const& old_linker : old_linkers)
  {
    auto old_linker_coo_groups = find_coo_groups(old_linker.atoms, old_linker.graph);
    auto new_transformed = kabsch_algorithm(new_linker.atoms, new_coo_groups, old_linker.atoms, old_linker_coo_groups);
    auto new_transformed_graph = calc_connect_graph_v2(new_transformed, new_transformed.size(), 1.3, "cordero");
    new_linkers.emplace_back(linker{.atoms=new_transformed, .graph=new_transformed_graph });
  }

  // HARDCODED OLD PARAMS
  double const a = 9.922710;
  double const b = 29.951300;
  double const c = 7.450510;
  double const alpha = 89.905197;
  double const beta = 132.923996;
  double const gamma = 90.234299;

  // Old linker len (same as between 1 and 3 atom in linker0)
  double const old_length = 9.91;
  double const delta = new_length - old_length;

  double const new_b = b + 2*delta;


  // GENERATING NEW CIF FINALLY
  // Making it like a sandwich

  cif_data new_cif;
  new_cif.cell = make_unit_cell(a, new_b, c, alpha, beta, gamma);

  // Layer 1 : lower node
  append_vector(old_lower_node, new_cif.abs_atoms);


  // Layer 2 : first linker
  double old_linker_lowest_point0 = min_y_coordinate(old_linkers[0].atoms);
  double new_linker_lowest_point0 = min_y_coordinate(new_linkers[0].atoms);
  double diff0 = new_linker_lowest_point0 - old_linker_lowest_point0;

  for(auto& atom : new_linkers[0].atoms) { atom.coord[1] += diff0 + delta; }
  auto new_linker0_noo = delete_o_from_coo(new_linkers[0].atoms, new_linkers[0].graph, new_coo_groups);
  append_vector(new_linker0_noo, new_cif.abs_atoms);

  // Layer 3: middle node
  for(auto& atom : old_middle_node) { atom.coord[1] += delta; }
  append_vector(old_middle_node, new_cif.abs_atoms);

  // Layer 4: third linker
  double old_linker_lowest_point1 = min_y_coordinate(old_linkers[1].atoms);
  double new_linker_lowest_point1 = min_y_coordinate(new_linkers[1].atoms);
  double diff1 = new_linker_lowest_point1 - old_linker_lowest_point1;

  for(auto& atom : new_linkers[1].atoms) { atom.coord[1] += diff1 + 2*delta; }
  auto new_linker1_noo = delete_o_from_coo(new_linkers[1].atoms, new_linkers[1].graph, new_coo_groups);
  append_vector(new_linker1_noo, new_cif.abs_atoms);

  // Layer 5: upper node
  for(auto& atom : old_upper_node) { atom.coord[1] += 2*delta; }
  append_vector(old_upper_node, new_cif.abs_atoms);

  std::ofstream final_xyz_ofs("lim.xyz");
  print_xyz(new_cif.abs_atoms, "qwe", &final_xyz_ofs);

  // POST PROCESSING

  // auto new_end = std::remove_if(std::begin(new_cif.abs_atoms), std::end(new_cif.abs_atoms),
  //                               [&new_cif](auto&& a){ return is_inside_cell(a, new_cif.cell); });
  // new_cif.abs_atoms.erase(new_end, std::end(new_cif.abs_atoms));
  for(auto& atom : new_cif.abs_atoms) { reflect_if_needed(atom, new_cif.cell); }
  new_cif.calc_fract_coords();
  std::ofstream final_cif_ofs("LIM.cif");
  print_cif(new_cif, &final_cif_ofs);

}
