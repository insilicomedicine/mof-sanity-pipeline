// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "mutation_ops.hpp"
#include "cif.hpp"
#include "cif_ops.hpp"
#include "data_tables.hpp"
#include "obabel_wrappers.hpp"
#include "rdkit_wrappers.hpp"
#include "utils.hpp"
#include "graph_ops.hpp"
#include "linalg.hpp"
#include <iterator>
#include <optional>

void remove_same_linkers(std::vector<linker>& linkers)
{
  auto linker_to_unique_indices = [](linker const& linker)
  {
      return std::set<std::size_t>(std::begin(linker.orig_indices), std::end(linker.orig_indices));
  };


  std::sort(std::begin(linkers), std::end(linkers),
            [&linker_to_unique_indices](linker const& lhs, linker const& rhs){ return linker_to_unique_indices(lhs) < linker_to_unique_indices(rhs); });
  auto new_end = std::unique(std::begin(linkers), std::end(linkers),
                             [&linker_to_unique_indices](linker const& lhs, linker const& rhs){ return linker_to_unique_indices(lhs) == linker_to_unique_indices(rhs); });
  linkers.erase(new_end, std::end(linkers));
}

std::vector<std::size_t> find_equivalent_atoms(std::size_t target_index, linker const& target_linker, std::vector<linker> const& group)
{
  std::vector<std::size_t> equiv_atoms;

  std::size_t canon_number = target_linker.canon_order[target_index];

  for(auto const& linker : group)
  {
    std::size_t needed_index = std::distance(std::begin(linker.canon_order),
                                             std::find(std::begin(linker.canon_order), std::end(linker.canon_order), canon_number));
    equiv_atoms.push_back(linker.orig_indices[needed_index]);
  }
  return equiv_atoms;
}

gene make_equivalent_atoms(gene const& old_gene, linker const& target_linker, std::vector<linker> const& group)
{
  if(old_gene.sites.size() != 1) { throw std::runtime_error("Wrong gene!!!!!!!"); }
  gene new_gene {old_gene.mutation_number, {}};
  std::vector<std::vector<std::size_t>> transpose;
  for(std::size_t old_site : old_gene.sites[0])
  {
    auto new_positions = find_equivalent_atoms(old_site, target_linker, group);
    transpose.push_back(new_positions);
  }

  new_gene.sites = std::vector<std::vector<std::size_t>>(transpose[0].size(), std::vector<std::size_t>(transpose.size()));
  for(std::size_t i=0; i<transpose.size(); ++i)
  {
    for(std::size_t j=0; j<transpose[0].size(); ++j)
    {
      new_gene.sites[j][i] = transpose[i][j];
    }
  }

  return new_gene;
}

std::optional<std::vector<std::size_t>> match_61(std::size_t target, linker const& linker)
{
  if(linker.atoms[target].name != "C") { return std::nullopt; }
  std::vector<std::size_t> ns_CH_group;
  int hydrogen = -1;
  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    if(linker.atoms[bond.neighbour].name == "H") { hydrogen = bond.neighbour; }
    else { ns_CH_group.push_back(bond.neighbour); }
  }
  if(hydrogen < 0) { return std::nullopt; }
  ns_CH_group.push_back(target);
  ns_CH_group.push_back(hydrogen);
  return ns_CH_group;
}

std::optional<std::vector<std::size_t>> match_62(std::size_t target, linker const& linker)
{
  if(linker.atoms[target].name != "C") { return std::nullopt; }
  std::vector<std::size_t> ns_CF_group;
  int fluorine = -1;
  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    if(linker.atoms[bond.neighbour].name == "F") { fluorine = bond.neighbour; }
    else { ns_CF_group.push_back(bond.neighbour); }
  }
  if(fluorine < 0) { return std::nullopt; }
  ns_CF_group.push_back(target);
  ns_CF_group.push_back(fluorine);
  return ns_CF_group;
}

std::optional<std::vector<std::size_t>> match_63(std::size_t target, linker const& linker)
{
  if(linker.atoms[target].name != "C") { return std::nullopt; }
  std::vector<std::size_t> ns_CCH3_group;
  std::vector<std::size_t> CH3_group {};
  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    if(is_CH3_group(bond.neighbour, linker.atoms, linker.graph))
    {
      CH3_group.push_back(bond.neighbour);
      for(auto const& target_bond : linker.graph.adjacency_list[bond.neighbour])
      {
        if(linker.atoms[target_bond.neighbour].name == "H")
        {
          CH3_group.push_back(target_bond.neighbour);
        }
      }
    }
    else { ns_CCH3_group.push_back(bond.neighbour); }
  }
  if(CH3_group.size() == 0) { return std::nullopt; }
  ns_CCH3_group.push_back(target);
  append_vector(CH3_group, ns_CCH3_group);
  return ns_CCH3_group;
}

std::optional<std::vector<std::size_t>> match_64(std::size_t target, linker const& linker)
{
  if(linker.atoms[target].name != "C") { return std::nullopt; }
  std::vector<std::size_t> ns_CNH2_group;
  std::vector<std::size_t> NH2_group {};
  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    if(is_NH2_group(bond.neighbour, linker.atoms, linker.graph))
    {
      NH2_group.push_back(bond.neighbour);
      for(auto const& target_bond : linker.graph.adjacency_list[bond.neighbour])
      {
        if(linker.atoms[target_bond.neighbour].name == "H")
        {
          NH2_group.push_back(target_bond.neighbour);
        }
      }
    }
    else { ns_CNH2_group.push_back(bond.neighbour); }
  }
  if(NH2_group.size() == 0) { return std::nullopt; }
  ns_CNH2_group.push_back(target);
  append_vector(NH2_group, ns_CNH2_group);
  return ns_CNH2_group;
}

std::optional<std::vector<std::size_t>> match_65(std::size_t target, linker const& linker)
{
  if(linker.atoms[target].name != "N") { return std::nullopt; }
  if(linker.graph.adjacency_list[target].size() != 2) { return std::nullopt; }
  std::vector<std::size_t> ns_N_group;
  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    ns_N_group.push_back(bond.neighbour);
  }
  ns_N_group.push_back(target);
  return ns_N_group;
}

gene match_atom_group_3(std::size_t target, linker const& linker)
{
  if(auto opt = match_61(target,linker); opt.has_value())
  {
    return gene{.mutation_number=61, .sites={opt.value()}};
  }
  if(auto opt = match_62(target,linker); opt.has_value())
  {
    return gene{.mutation_number=62, .sites={opt.value()}};
  }
  if(auto opt = match_63(target,linker); opt.has_value())
  {
    return gene{.mutation_number=63, .sites={opt.value()}};
  }
  if(auto opt = match_64(target,linker); opt.has_value())
  {
    return gene{.mutation_number=64, .sites={opt.value()}};
  }
  if(auto opt = match_65(target,linker); opt.has_value())
  {
    return gene{.mutation_number=65, .sites={opt.value()}};
  }
  return gene{-1, {}};
}

std::vector<gene> find_sites_group_1(std::vector<atom> const& atoms)
{
  std::vector<gene> genes;
  for(auto const& [metal, num] : data::mutations::group_1)
  {
    gene curr_gene {.mutation_number=num, .sites{}};
    std::vector<std::size_t> sites;
    for(std::size_t i=0; i<atoms.size(); ++i)
    {
      if(atoms[i].name == metal) { sites.push_back(i); }
    }
    if(sites.size() != 0) { curr_gene.sites.push_back(sites); }
    if(curr_gene.sites.size() != 0) { genes.push_back(curr_gene); }
  }
  return genes;
}

std::vector<gene> find_sites_group_2(std::vector<atom> const& atoms)
{
  std::vector<gene> genes;
  for(auto const& [metal, num] : data::mutations::group_2)
  {
    gene curr_gene {.mutation_number=num, .sites{}};
    std::vector<std::size_t> sites;
    for(std::size_t i=0; i<atoms.size(); ++i)
    {
      if(atoms[i].name == metal) { sites.push_back(i); }
    }
    if(sites.size() != 0) { curr_gene.sites.push_back(sites); }
    if(curr_gene.sites.size() != 0) { genes.push_back(curr_gene); }
  }
  return genes;
}

std::vector<gene> find_sites_group_3(cif_data& cif,
                                     double r_min, double r_max,
                                     double tolerance,
                                     std::string const& strategy)
{
  auto linkers = separate_linkers(cif, r_min, r_max, tolerance, strategy);
  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.smiles = obabel_wrappers::atoms_to_canon_smiles(linker.atoms);
    linker.graph = calc_connect_graph(linker.atoms, linker.atoms.size(), tolerance, strategy);
    linker.orig_indices = transform_to_original_indices(linker.indices, cif.num_atoms);
    linker.canon_order = rdkit_wrappers::calc_canon_order_rdkit(linker.atoms);
  }

  remove_same_linkers(linkers);
  auto linker_groups = group_by_smiles(linkers);

  std::vector<gene> genes;

  for(auto const& [smiles, group] : linker_groups)
  {
    auto rings = find_all_n_member_rings(group[0].graph, 6);
    for(auto const& ring : rings)
    {
      if(!all_atoms_unsaturated(ring, group[0].atoms, group[0].graph)) { continue; }
      for(std::size_t index : ring)
      {
        auto tmp_gene = match_atom_group_3(index, group[0]);
        if(tmp_gene.mutation_number > 0)
        {
          genes.push_back(make_equivalent_atoms(tmp_gene, group[0], group));
        }
      }
    }
  }
  return genes;
}

bool is_group_1(int code) { return (code>0) and (code < 7); }
bool is_group_2(int code) { return (code>20) and (code < 25); }
bool is_group_3(int code) { return (code>60) and (code < 66); }

void mutate_group_1(std::size_t target, int new_mutation, cif_data& cif)
{
  if(!is_group_1(cif.cmsm.genome[target].mutation_number)) { throw std::domain_error("Wrong mutation provided!!"); }
  if(new_mutation == cif.cmsm.genome[target].mutation_number) { return; }
  for([[maybe_unused]] auto site : cif.cmsm.genome[target].sites)
  {
  }
}

void mutate_group_2(std::size_t target, int new_mutation, cif_data& cif)
{
  if(!is_group_2(cif.cmsm.genome[target].mutation_number)) { throw std::domain_error("Wrong mutation provided!!"); }
  if(new_mutation == cif.cmsm.genome[target].mutation_number) { return; }
  for([[maybe_unused]] auto site: cif.cmsm.genome[target].sites)
  {
  }
}

atom make_atom_aroma(atom const& neighbour1, atom const& neighbour2, atom const& center, std::string const& name)
{
  vec3d n1c_vec = vec3d(center.coord) - vec3d(neighbour1.coord);
  vec3d n2c_vec = vec3d(center.coord) - vec3d(neighbour2.coord);

  vec3d direction_vec = n1c_vec + n2c_vec;
  normalize(direction_vec);

  double new_length = data::cov_radii_table["C"].covalent_radius_cordero/100.0 + data::cov_radii_table[name].covalent_radius_cordero/100.0;

  vec3d new_atom_coord = vec3d(center.coord) + direction_vec*new_length;
  atom new_atom {.coord=new_atom_coord.to_array(), .name=name};

  return new_atom;
}

std::vector<atom> make_CH3_aroma(atom const& neighbour1, atom const& neighbour2, atom const& center)
{
  vec3d n1c_vec = vec3d(center.coord) - vec3d(neighbour1.coord);
  vec3d n2c_vec = vec3d(center.coord) - vec3d(neighbour2.coord);

  vec3d direction_vec = n1c_vec + n2c_vec;
  normalize(direction_vec);

  double cc_length = data::cov_radii_table["C"].covalent_radius_cordero*2.0/100.0;

  vec3d new_c_coord = vec3d(center.coord) + direction_vec*cc_length;
  atom carbon {.coord=new_c_coord.to_array(), .name="C"};

  double const len_to_plane = 0.62735682157; // 1.087 * cos(54.75deg)
  double const len_to_h = 0.887689370461; // 1.087 * sin(54.75deg)

  vec3d plane_start_vec = new_c_coord + direction_vec*len_to_plane;
  vec3d first_h = {-direction_vec[1], direction_vec[0], 0};
  vec3d second_h = rotate_around_axis(first_h, direction_vec, 120*M_PI/180);
  vec3d third_h = rotate_around_axis(second_h, direction_vec, 120*M_PI/180);

  first_h = first_h*len_to_h + plane_start_vec;
  second_h = second_h*len_to_h + plane_start_vec;
  third_h = third_h*len_to_h + plane_start_vec;

  atom h1 {.coord=first_h.to_array(), .name="H"};
  atom h2 {.coord=second_h.to_array(), .name="H"};
  atom h3 {.coord=third_h.to_array(), .name="H"};

  return {carbon, h1, h2, h3};
}

void mutate_61(std::size_t target, cif_data& cif)
{
  for(auto const& sites : cif.cmsm.genome[target].sites)
  {
    if(sites.size() < 4)
    {
      atom new_hydrogen = make_atom_aroma(cif.abs_atoms[sites[0]], cif.abs_atoms[sites[1]], cif.abs_atoms[sites[2]], "H");
      reflect_if_needed(new_hydrogen, cif);
      cif.abs_atoms.push_back(new_hydrogen);
      continue;
    }
    if(sites.size() > 4)
    {
      std::cerr << "ERASING ATOMS NOT NEAR C\n";
      for(std::size_t i=4; i<sites.size(); ++i)
      {
        cif.abs_atoms.erase(std::next(std::begin(cif.abs_atoms), sites[i]));
      }
    }
    cif.abs_atoms[sites[3]].name = "H";
  }
}

void mutate_62(std::size_t target, cif_data& cif)
{
  for(auto const& sites : cif.cmsm.genome[target].sites)
  {
    if(sites.size() < 4)
    {
      atom new_fluorine = make_atom_aroma(cif.abs_atoms[sites[0]], cif.abs_atoms[sites[1]], cif.abs_atoms[sites[2]], "F");
      reflect_if_needed(new_fluorine, cif);
      cif.abs_atoms.push_back(new_fluorine);
      continue;
    }
    if(sites.size() > 4)
    {
      std::cerr << "ERASING ATOMS NOT NEAR C\n";
      for(std::size_t i=4; i<sites.size(); ++i)
      {
        cif.abs_atoms.erase(std::next(std::begin(cif.abs_atoms), sites[i]));
      }
    }
    cif.abs_atoms[sites[3]].name = "F";
  }
}

void mutate_63(std::size_t target, cif_data& cif)
{
  for(auto const& sites : cif.cmsm.genome[target].sites)
  {
    if(sites.size() > 3)
    {
      std::cerr << "ERASING ATOMS NOT NEAR C\n";
      for(std::size_t i=3; i<sites.size(); ++i)
      {
        cif.abs_atoms.erase(std::next(std::begin(cif.abs_atoms), sites[i]));
      }
    }
    cif.abs_atoms[sites[2]].name = "C";
    std::vector<atom> new_CH3_group = make_CH3_aroma(cif.abs_atoms[sites[0]], cif.abs_atoms[sites[1]], cif.abs_atoms[sites[2]]);
    for(auto& atom : new_CH3_group) { reflect_if_needed(atom, cif); }
    append_vector(new_CH3_group, cif.abs_atoms);
  }
}

void mutate_64([[maybe_unused]] std::size_t target, [[maybe_unused]] cif_data& cif)
{
  throw std::runtime_error("not implemented");
}

void mutate_65(std::size_t target, cif_data& cif)
{
  for(auto const& sites : cif.cmsm.genome[target].sites)
  {
    if(sites.size() > 3)
    {
      std::cerr << "ERASING ATOMS EXCEPT N\n";
      for(std::size_t i=3; i<sites.size(); ++i)
      {
        cif.abs_atoms.erase(std::next(std::begin(cif.abs_atoms), sites[i]));
      }
    }
    cif.abs_atoms[sites[2]].name = "N";
  }
}

void mutate_group_3(std::size_t target, int new_mutation, cif_data& cif)
{
  if(!is_group_3(cif.cmsm.genome[target].mutation_number)) { throw std::domain_error("Wrong mutation provided!!"); }
  if(new_mutation == cif.cmsm.genome[target].mutation_number) { return; }

  switch(new_mutation)
  {
    case 61:
      mutate_61(target, cif);
      break;
    case 62:
      mutate_62(target, cif);
      break;
    case 63:
      mutate_63(target, cif);
      break;
    case 64:
      mutate_64(target, cif);
      break;
    case 65:
      mutate_65(target, cif);
      break;
  }
}
