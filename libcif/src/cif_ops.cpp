// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif_ops.hpp"
#include "cif.hpp"
#include "data_tables.hpp"
#include "utils.hpp"
#include "graph_ops.hpp"

#include <iterator>
#include <regex>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <map>
#include <unordered_set>

double calc_surface(std::vector<atom_with_radius> const& spheres)
{
	double S = 0.0;
	for(auto const& sphere : spheres)
	{
    S += sphere_surface(sphere.vdw_radius);
	}
	return S;
}

std::vector<atom_with_radius> select_polar_atoms(std::vector<atom_with_radius> const& spheres, double thr)
{
	std::vector<atom_with_radius> polar;
	std::copy_if(std::begin(spheres), std::end(spheres),
				 std::back_inserter(polar),
				 [thr](auto&& sphere){ return std::abs(sphere.q) > thr; }
				 );
	return polar;
}

std::vector<std::size_t> find_metal_idxs(std::vector<atom> const& atoms)
{
	std::vector<std::size_t> metal_idxs;
	for(std::size_t i=0; i<atoms.size(); ++i)
	{
		if(data::is_metal(atoms[i])) { metal_idxs.push_back(i); }
	}
	return metal_idxs;
}

std::vector<std::size_t> find_atoms_closer_than(std::size_t target, std::vector<atom> const& atoms, double cutoff)
{
	std::vector<std::size_t> res;
	for(std::size_t i=0; i<atoms.size(); ++i)
	{
		if(i==target) { continue; }
		
		double const r = atoms_range(atoms[target], atoms[i]);
		if(r < cutoff) { res.push_back(i); }
	}
	return res;
}

double closest_neighbour_range(std::size_t target, std::vector<atom> const& atoms)
{
	std::vector<double> ranges(atoms.size());
	for(std::size_t i=0; i<atoms.size(); ++i)
	{
		ranges[i] = atoms_range(atoms[target], atoms[i]);
	}
	double closest_range = std::numeric_limits<double>::max();
	for(std::size_t i=0; i<atoms.size(); ++i)
	{
		if(i==target) { continue; }
		if(ranges[i] < closest_range) { closest_range = ranges[i]; }
	}
	return closest_range;
}

// If only metal neighbours, than return max valence,
// otherwise return valence as if all metals are ONE atom
std::size_t reduce_metal_neighbours(std::size_t target, std::size_t curr_valence, cif_data const& cif)
{
  std::size_t num_metals = std::count_if(std::begin(cif.graph.adjacency_list[target]), std::end(cif.graph.adjacency_list[target]),
                                          [&cif](auto&& bond){ return data::is_metal(cif.supercell[bond.neighbour]); });
  if(num_metals == 0) { return curr_valence; }
  std::cerr << "BRIDGE ATOM " << cif.abs_atoms[target].name << " #" << target+1 << " has " << num_metals << " metal neighbours" << std::endl;
  if(curr_valence == num_metals) { return data::get_max_valence(cif.abs_atoms[target].name); }
  else { return curr_valence - num_metals; }
}

bool cif_is_good(cif_data& cif,
                 double r_min, double r_max,
                 double tolerance,
                 std::string const& strategy)
{
  if(!calc_graph(cif, r_min, r_max, tolerance, strategy)) { return false; }
  if(!check_disconnected(cif)) { return false; }

  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    if(data::is_metal(cif.abs_atoms[i])) { continue; }

    std::size_t num_valence = cif.graph.adjacency_list[i].size();
    std::size_t max_valence = data::get_max_valence(cif.abs_atoms[i].name);

    if(num_valence > max_valence)
    {
      std::size_t new_valence = num_valence;
      if(data::is_bridge_atom(cif.abs_atoms[i].name))
      {
        new_valence = reduce_metal_neighbours(i, num_valence, cif);
      }
      if(new_valence > max_valence)
      {
        std::cerr << "Atom " << cif.abs_atoms[i].name << " #" << i+1 << " has more than " << max_valence << " neighbours!" << std::endl;
        return false;
      }
    }
  }

	return true;
}

bool cif_is_good_v2(cif_data& cif,
                    double r_min, double r_max,
                    double tolerance,
                    std::string const& strategy)
{
  if(!calc_graph(cif, r_min, r_max, tolerance, strategy)) { return false; }
  if(!check_disconnected(cif)) { return false; }

  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    if(data::is_metal(cif.abs_atoms[i])) { continue; }

    std::size_t num_valence = cif.graph.adjacency_list[i].size();
    std::size_t max_valence = data::get_max_valence(cif.abs_atoms[i].name);

    if(num_valence > max_valence)
    {
      std::size_t new_valence = num_valence;
      if(data::is_bridge_atom(cif.abs_atoms[i].name))
      {
        new_valence = reduce_metal_neighbours(i, num_valence, cif);
      }
      if(new_valence > max_valence)
      {
        std::cerr << "Atom " << cif.abs_atoms[i].name << " #" << i+1 << " has more than " << max_valence << " neighbours!" << std::endl;
        return false;
      }
    }
  }

  if(!check_carbons(cif.graph, cif.supercell))
  {
    std::cerr << "WARNING! Cif **CARBON** angles are not optimal!" << std::endl;
  }
  if(!check_nitrogens(cif.graph, cif.supercell))
  {
    std::cerr << "WARNING! Cif **NITROGEN** angles are not optimal!" << std::endl;
  }

	return true;
}

double cif_quality(cif_data& cif,
                 double r_min, double r_max,
                 double tolerance,
                 std::string const& strategy,
                 double penalty)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1);
  cif.graph = calc_connect_graph(cif.supercell, cif.num_atoms, tolerance, strategy);

  // Quality 1 -- number of overlapping atoms
  auto overlaps = find_overlap_atoms(cif.graph, r_min);

  remove_too_long_bonds(cif.graph, r_max);

  // Quality 2 -- number of too many neighbours or too few
  std::vector<std::size_t> too_many;
  std::vector<std::size_t> special;
  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    if(data::is_metal(cif.abs_atoms[i])) { continue; }

    std::size_t num_valence = cif.graph.adjacency_list[i].size();
    std::size_t max_valence = data::get_max_valence(cif.abs_atoms[i].name);

    if(num_valence > max_valence)
    {
      std::size_t new_valence = num_valence;
      if(data::is_bridge_atom(cif.abs_atoms[i].name))
      {
        new_valence = reduce_metal_neighbours(i, num_valence, cif);
      }
      if(new_valence > max_valence)
      {
        too_many.push_back(i);
      }
    }
    else if(num_valence < max_valence)
    {
      if(std::all_of(std::begin(cif.graph.adjacency_list[i]), std::end(cif.graph.adjacency_list[i]),
                     [&](auto&& bond){ return !data::is_unsaturated(cif.supercell[bond.neighbour].name, num_valence); }))
      {
        special.push_back(i);
      }
    }
  }

  std::unordered_set<std::size_t> unique_bad_atoms;
  for(auto idx : overlaps) { unique_bad_atoms.insert(idx); }
  for(auto idx : too_many) { unique_bad_atoms.insert(idx); }
  for(auto idx : special) { unique_bad_atoms.insert(idx); }

  double uniq_bad = unique_bad_atoms.size();
  double all = cif.abs_atoms.size();

  return std::pow(uniq_bad/all, penalty);
}

void make_linker_dfs(std::size_t target, cif_data const& cif,
                     std::vector<std::size_t>& linker,
                     std::unordered_set<std::size_t>& visited,
                     double tolerance, std::string const& strategy)
{
  if(data::is_metal(cif.supercell[target])) { return; }
	linker.push_back(target);
	visited.insert(target);

  auto const& bonds = find_bonds(target, cif.supercell, tolerance, strategy);
  for(auto const& bond : bonds)
  {
    if(visited.contains(bond.neighbour)) { continue; }
    make_linker_dfs(bond.neighbour, cif, linker, visited, tolerance, strategy);
  }
}

void make_node_dfs(std::size_t target, cif_data const& cif,
                   std::vector<std::size_t>& node,
                   std::unordered_set<std::size_t>& visited,
                   std::unordered_set<std::size_t> const& linkers_indices,
                   double tolerance, std::string const& strategy)
{
  if(linkers_indices.contains(target)) { return; }
  node.push_back(target);
  visited.insert(target);

  auto const& bonds = find_bonds(target, cif.supercell, tolerance, strategy);
  for(auto const& bond : bonds)
  {
    if(visited.contains(bond.neighbour)) { continue; }
    make_node_dfs(bond.neighbour, cif, node, visited, linkers_indices, tolerance, strategy);
  }
}

std::vector<atom> indices_to_atoms(std::vector<std::size_t> const& indices, std::vector<atom> const& orig_atoms)
{
  std::vector<atom> atoms;
  atoms.reserve(indices.size());
  for(auto index : indices)
  {
    atoms.push_back(orig_atoms[index]);
  }
  return atoms;
}




/*
 * Now this function does the following:
 *
 * 1) Finds all metal atoms and all its neighbours
 * 2) Calls linker finding algorithm from EACH neighbour
 *  for each metal
 *
 * Comments:
 * - The algorithm stores already visited atoms, so
 *   each atoms can only belong to a single linker
 * - The algorithm is started only for atoms in the
 *   original cell (not generated 26)
 */
std::vector<linker> separate_linkers(cif_data& cif,
                                     double r_min, double r_max,
                                     double tolerance,
                                     std::string const& strategy)
{
  if(!calc_graph(cif, r_min, r_max, tolerance, strategy)) { throw std::runtime_error("CIF HAS MISTAKES!!"); }
  if(!check_disconnected(cif)) { std::cerr << "WARNING!! Cif contains disconnected subgraphs!"; }

	std::vector<linker> linkers;
	std::unordered_set<std::size_t> visited;


	auto metal_idxs = find_metal_idxs(cif.abs_atoms);
	if(metal_idxs.size() == 0) { throw std::runtime_error("No metals found! Can't separate linker!"); }

	for(std::size_t metal_idx : metal_idxs)
	{
    auto linker_starts = find_bonds(metal_idx, cif.supercell, tolerance, strategy);
    for(auto const& bond : linker_starts)
    {
      if(visited.contains(bond.neighbour)) { continue; }
      if(data::is_metal(cif.supercell[bond.neighbour])) { continue; }
      std::vector<std::size_t> linker_indices;
      make_linker_dfs(bond.neighbour, cif, linker_indices, visited, tolerance, strategy);
      if(linker_indices.size() < 2) { std::cerr << "WARNING! Found linker with size less than 2!" << std::endl; continue; }
      linkers.push_back({.indices=linker_indices});
    }
	}

  if(linkers.size() == 0) { throw std::runtime_error("ERROR, WRITE ME TO FIX IT!!!"); }
  return linkers;
}

std::vector<linker> separate_linkers_v2(cif_data &cif, std::size_t min_size,
                                        double r_min, double r_max,
                                        double tolerance,
                                        std::string const& strategy)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1);
  cif.graph = calc_connect_graph_v2(cif.supercell, cif.supercell.size(), tolerance, strategy);
  if(!verify_trim_graph(cif.graph, cif.supercell, r_min, r_max)) { throw std::runtime_error("Graph is BAD!"); }

  auto metal_idxs = find_metal_idxs(cif.abs_atoms);
  if(metal_idxs.size() == 0) { throw std::runtime_error("No metals found!"); }

  std::vector<linker> linkers;
  std::unordered_set<std::size_t> visited;

  for(std::size_t metal_idx : metal_idxs)
  {
    for(auto const& bond : cif.graph.adjacency_list[metal_idx])
    {
      if(visited.contains(bond.neighbour)) { continue; }
      if(data::is_metal(cif.supercell[bond.neighbour])) { continue; }
      std::vector<std::size_t> linker_indices;
      make_linker_dfs_v2(bond.neighbour, cif.graph, cif.supercell, linker_indices, visited);
      if(linker_indices.size() < min_size) { std::cerr << "WARNING! Found linker with size less than " << min_size << std::endl; continue;}
      linkers.push_back({.indices=linker_indices});
    }
  }

  if(linkers.size() == 0) { throw std::runtime_error("ERROR, NO LINKERS FOUND!"); }
  return linkers;
}

bool check_infinite(std::vector<std::size_t> const& indices, std::size_t orig_size)
{
  auto orig_indices = transform_to_original_indices(indices, orig_size);
  std::unordered_set<std::size_t> unique_indices;
  for(std::size_t index : orig_indices)
  {
    if(unique_indices.contains(index)) { return true; }
    unique_indices.insert(index);
  }
  return false;
}

bool has_no_orig_indices(std::vector<std::size_t> indices, std::size_t orig_size)
{
  auto is_orig = [orig_size](std::size_t idx){ return idx < orig_size; };
  return !std::any_of(std::begin(indices),
                      std::end(indices),
                      is_orig);
}

void remove_all_outside(std::vector<linker>& structs, std::size_t orig_size)
{
  auto new_end = std::remove_if(std::begin(structs),
                                std::end(structs),
                                [orig_size](auto&& l)
                                { return has_no_orig_indices(l.indices, orig_size); });
  structs.erase(new_end, std::end(structs));
}

bool is_me_pillar(std::size_t target, molecular_graph const& graph, std::vector<atom> const& atoms)
{
  if(!data::is_pillar_me(atoms[target])) { return false; }

  if(graph.adjacency_list[target].size() != 5 and
     graph.adjacency_list[target].size() != 6) { return false; }

  auto is_f_or_o = [&atoms](bond const& b)
  {
    return atoms[b.neighbour].name == "F" or atoms[b.neighbour].name == "O";
  };

  if(!std::ranges::all_of(graph.adjacency_list[target], is_f_or_o)) { return false; }
  return std::find_if(std::begin(graph.adjacency_list[target]), std::end(graph.adjacency_list[target]),
                      [&atoms](bond const& b){ return atoms[b.neighbour].name == "F"; }) != std::end(graph.adjacency_list[target]);
}

std::string calc_cif_formula(cif_data const& cif)
{
	std::vector<std::string> atom_names;
  auto const& some_atoms = cif.get_any_atoms_ref();
	atom_names.reserve(some_atoms.size());
	std::transform(std::begin(some_atoms), std::end(some_atoms),
			std::back_inserter(atom_names), [](atom const& a){ return a.name; });
  std::sort(std::begin(atom_names), std::end(atom_names));

	std::vector<std::string> unique_atom_names(atom_names);
	auto last_it = std::unique(std::begin(unique_atom_names), std::end(unique_atom_names));
	unique_atom_names.erase(last_it, std::end(unique_atom_names));

	std::string formula{};
	for(std::string const& name : unique_atom_names)
	{
		std::size_t num = std::count(std::begin(atom_names), std::end(atom_names), name);
		formula += name;
    formula += '_';
		formula += std::to_string(num);
    formula += '_';
	}
	formula.pop_back();
	return formula;
}

double calc_volume(cif_data const& cif)
{
  double const a = cif.cell.a;
  double const b = cif.cell.b;
  double const c = cif.cell.c;
  double const cos_alpha = std::cos(ang_to_rad(cif.cell.alpha));
  double const cos_beta = std::cos(ang_to_rad(cif.cell.beta));
  double const cos_gamma = std::cos(ang_to_rad(cif.cell.gamma));

  double V = 1 + 2*cos_alpha*cos_beta*cos_gamma;
  V -= cos_alpha*cos_alpha + cos_beta*cos_beta + cos_gamma*cos_gamma;
  V = std::sqrt(V);
  V *= a*b*c;

  return V;
}

double calc_mass(cif_data const& cif)
{
  double m {0};
  auto const& atoms = cif.get_any_atoms_ref();
  for(auto const& atom : atoms)
  {
    m += data::get_mass(atom.name);
  }
  return m;
}

// returns density in g/cm^3
double calc_density(cif_data const& cif)
{
  double const mass = calc_mass(cif) * amu_to_g;
  double const volume = calc_volume(cif) * angstrom3_to_cm3;
  return mass/volume;
}

double calc_PSA_ratio(cif_data const& cif, double thr, std::string const& method)
{
	std::vector<atom_with_radius> spheres;
	spheres.reserve(cif.num_atoms);
  auto const& some_atoms = cif.get_any_atoms_ref();
	for(atom const& a : some_atoms)
	{
		double vdw_radius = data::get_vdw_radius(a.name, method);
		spheres.emplace_back(a, vdw_radius);
	}
	double total_surface = calc_surface(spheres);


	auto polar_spheres = select_polar_atoms(spheres, thr);
	double polar_surface = calc_surface(polar_spheres);


	return polar_surface / total_surface;
}

double calc_total_charge(cif_data const& cif)
{
  auto const& atoms = cif.get_any_atoms_ref();
  double charge = 0;
  for(auto const& atom : atoms)
  {
    charge += atom.q;
  }
  return charge;
}

#ifdef LIBCIF_BUILD_OBABEL

#include "obabel_wrappers.hpp"

std::unordered_set<std::string> separate_linkers_smiles(cif_data& cif,
                                                double r_min, double r_max,
                                                double tolerance, std::string const& strategy)
{
  auto linkers = separate_linkers(cif, r_min, r_max, tolerance, strategy);
  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.smiles = obabel_wrappers::atoms_to_canon_smiles(linker.atoms);
  }
  std::unordered_set<std::string> smiless;
  for(auto const& linker : linkers)
  {
    smiless.insert(linker.smiles);
  }
  return smiless;
}

std::unordered_set<std::string> separate_linkers_smiles_v2(cif_data& cif, std::size_t min_size,
                                                           double r_min, double r_max,
                                                           double tolerance, std::string const& strategy)
{
  auto linkers = separate_linkers_v2(cif, min_size, r_min, r_max, tolerance, strategy);
  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.smiles = obabel_wrappers::atoms_to_canon_smiles(linker.atoms);
  }
  std::unordered_set<std::string> smiless;
  for(auto const& linker : linkers)
  {
    smiless.insert(linker.smiles);
  }
  return smiless;
}
#endif

#ifdef LIBCIF_BUILD_RDKIT

#include "rdkit_wrappers.hpp"
#include "mutation_ops.hpp"

chromosome make_chromosome(cif_data& cif,
                           double r_min, double r_max,
                           double tolerance,
                           std::string const& strategy)
{
  cif.calc_abs_coords();

  chromosome cmsm;
  append_vector(find_sites_group_1(cif.abs_atoms), cmsm.genome);
  append_vector(find_sites_group_2(cif.abs_atoms), cmsm.genome);
  append_vector(find_sites_group_3(cif, r_min, r_max, tolerance, strategy), cmsm.genome);

  return cmsm;
}

void print_mutant(cif_data const& cif, std::ostream* os)
{
  cif.dump_primcell(os);
  cif.cmsm.print_full(os);
}

cif_data read_mutant(std::istream& is)
{
  cif_data cif;
  // READ PRIMCELL PARAMS
  double a, b, c, alpha, beta, gamma;
  is >> a >> b >> c >> alpha >> beta >> gamma;
  cif.cell = make_unit_cell(a, b, c, alpha, beta, gamma);
  // READ XYZ ABS COORDS
  is >> cif.num_atoms;
  std::string tmp;
  std::getline(is, tmp);
  std::getline(is, tmp);
  cif.abs_atoms.resize(cif.num_atoms);
  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    is >> cif.abs_atoms[i].name >> cif.abs_atoms[i].coord[0] >> cif.abs_atoms[i].coord[1] >> cif.abs_atoms[i].coord[2];
  }

  chromosome cmsm;

  // READ FIRST LINE (OF CODES ONLY)
  std::vector<int> mutation_codes;
  std::getline(is, tmp);
  std::getline(is, tmp);
  std::istringstream iss(tmp);
  std::string tmp2;
  while(iss >> tmp2)
  {
    mutation_codes.push_back(std::stoi(tmp2));
  }
  cmsm.genome.resize(mutation_codes.size());
  
  for(std::size_t i=0; i<mutation_codes.size(); ++i)
  {
    int curr_mut_code;
    std::size_t num_sites;
    is >> curr_mut_code >> num_sites;
    if(mutation_codes[i] != curr_mut_code) { throw std::runtime_error("Wrong mutation code!"); }
    cmsm.genome[i].mutation_number = curr_mut_code;
    cmsm.genome[i].sites.resize(num_sites);
    std::getline(is, tmp2); // read the EOL
    for(std::size_t j=0; j<num_sites; ++j)
    {
      std::string site_atoms;
      std::getline(is, site_atoms);
      std::istringstream iss2(site_atoms);
      std::string single_atom;
      std::vector<std::size_t> site_indices;
      while(iss2 >> single_atom)
      {
        site_indices.push_back(std::stoi(single_atom));
      }
      cmsm.genome[i].sites[j] = site_indices;
    }
  }
  cif.cmsm = cmsm;

  return cif;
}

void mutate(cif_data& cif, std::vector<int> const& new_mutations)
{
  for(std::size_t i=0; i<new_mutations.size(); ++i)
  {
    if(is_group_1(new_mutations[i])) { mutate_group_1(i, new_mutations[i], cif); }
    else if(is_group_2(new_mutations[i])) { mutate_group_2(i, new_mutations[i], cif); }
    else if(is_group_3(new_mutations[i])) { mutate_group_3(i, new_mutations[i], cif); }
    else
    {
      throw std::domain_error(std::string{"ERROR!! Mutation doesn't belong to any group! num: "} + std::to_string(new_mutations[i]));
    }
  }
}

std::vector<double> determine_possible_charges(linker const& linker)
{
  if(linker.atoms.size() == 2)
  {
    // Check for free OH(1-)
    if((linker.atoms[0].name == "O" or linker.atoms[1].name == "O") and
       (linker.atoms[1].name == "H" or linker.atoms[0].name == "H"))
    {
      return {-1};
    }
    // Check for free NN(2-)
    else if(linker.atoms[0].name == "N" and linker.atoms[1].name == "N")
    {
      return {-2};
    }
    // Check for free CN group
    else if((linker.atoms[0].name == "C" or linker.atoms[1].name == "C") and
            (linker.atoms[1].name == "N" or linker.atoms[0].name == "N"))
    {
      return {-1};
    }
  }

  std::vector<double> possible_charges {};
  double curr_charge = 0;
  // Find COO group
  for(std::size_t i=0; i<linker.atoms.size(); ++i)
  {
    if(is_COO_group(i, linker.atoms, linker.graph))
    {
      curr_charge -= 1;
    }
    if(is_N_charged_group(i, linker.atoms, linker.graph))
    {
      curr_charge += 1;
    }
  }

  while(possible_charges.size() < 2)
  {
    try
    {
      rdkit_wrappers::determine_bonds(linker.atoms, curr_charge);
    }
    catch(std::exception& e)
    {
      curr_charge -= 1;
      continue;
    }
    std::cout << "Could be " << curr_charge << std::endl;
    possible_charges.push_back(curr_charge);
    curr_charge -= 1;
  }

  return possible_charges;
}

std::tuple<std::vector<std::string>,int> calc_metal_charges(cif_data& cif,
                                                            double r_min, double r_max,
                                                            double tolerance,
                                                            std::string const& strategy)
{
  double charge = 0;
  std::vector<std::string> metals {};

  auto linkers = separate_linkers_v2(cif, 0, r_min, r_max, tolerance, strategy);
  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.graph = calc_connect_graph(linker.atoms, linker.atoms.size(), tolerance, strategy);
    linker.smiles = rdkit_wrappers::atoms_to_canon_smiles(linker.atoms);
  }

  std::unordered_map<std::string, std::vector<double>> smiles_to_charges;
  for(auto const& linker : linkers)
  {
    if(smiles_to_charges.contains(linker.smiles)) { continue; }
    smiles_to_charges[linker.smiles] = determine_possible_charges(linker);
    for(auto charge : smiles_to_charges[linker.smiles]) { std::cout << charge << std::endl;}
  }

  for(auto& linker : linkers)
  {
    std::vector<double> possible_charges = determine_possible_charges(linker);
    if(possible_charges.size() == 1) { charge += possible_charges[0]; continue; }

    std::vector<std::size_t> metal_indexes;

    for(std::size_t index : linker.indices)
    {
      auto metal_it = std::find_if(std::begin(cif.graph.adjacency_list[index]),
                     std::end(cif.graph.adjacency_list[index]),
                      [](auto&& bond){ return bond.type == bond_type::METAL; });
      if(metal_it != std::end(cif.graph.adjacency_list[index]))
      {
        metal_indexes.push_back(metal_it->neighbour);
      }
    }

    std::size_t num_metals = metal_indexes.size();
    std::size_t num_outside = 0;

    for(std::size_t metal_index : metal_indexes)
    {
      if(metal_index >= cif.num_atoms) { num_outside += 1; }
    }
  }

  for(std::size_t i=0; i<cif.abs_atoms.size(); ++i)
  {
    if(data::is_metal(cif.supercell[i])) { metals.push_back(cif.supercell[i].name); }
    if(data::is_bridge_atom(cif.supercell[i].name))
    {
      if(std::all_of(cif.graph.adjacency_list[i].begin(),
                     cif.graph.adjacency_list[i].end(),
                     [&cif](auto&& bond){ return data::is_metal(cif.supercell[bond.neighbour]); }))
      {
        charge += data::bridge_charges[cif.supercell[i].name];
      }
    }
    if(data::is_pillar_cation(cif.supercell[i].name))
    {
      if(cif.graph.adjacency_list[i].size() == 6)
      {
        metals.pop_back();
        charge += data::pillar_cations[cif.supercell[i].name];
      }
    }
  }


  return {metals, charge};
}

#endif

#ifdef LIBCIF_BUILD_SPGLIB

#include "spglib_wrappers.hpp"

void cell_prim_spg(cif_data& cif, double tolerance)
{
	spglib_wrappers::primitivize_cell(cif, tolerance);
}

#endif

#ifdef LIBCIF_BUILD_HASH

void remove_subsets(std::vector<linker>& linkers)
{
  // Sort by size
  std::ranges::sort(linkers,
                    [](linker const& l1, linker const& l2){ return l1.atoms.size() < l2.atoms.size(); });
  // Transform to sets
  std::vector<std::unordered_set<std::size_t>> orig_indices_sets(linkers.size());
  std::transform(std::begin(linkers), std::end(linkers), std::begin(orig_indices_sets),
                 [](linker const& l)
                 { return std::unordered_set<std::size_t>(std::begin(l.orig_indices), std::end(l.orig_indices)); });
  for(std::size_t i=0; i<linkers.size(); ++i)
  {
    for(std::size_t j=i+1; j<linkers.size(); ++j)
    {
      if(is_subset(orig_indices_sets[i], orig_indices_sets[j]))
      {
        linkers[i].atoms.clear();
      }
    }
  }

  auto new_end = std::remove_if(std::begin(linkers), std::end(linkers),
                                [](linker const& l){ return l.atoms.size() == 0; });
  linkers.erase(new_end, std::end(linkers));
}

void remove_outside_dfs(std::size_t target,
                        linker const& linker,
                        std::vector<std::size_t>& keep_indices,
                        std::unordered_set<std::size_t>& visited,
                        std::unordered_set<std::size_t>& orig_visited)
{
  visited.insert(target);
  orig_visited.insert(linker.orig_indices[target]);
  keep_indices.push_back(target);

  for(auto const& bond : linker.graph.adjacency_list[target])
  {
    if(visited.contains(bond.neighbour)) { continue; }
    if(orig_visited.contains(linker.orig_indices[bond.neighbour])) { continue; }
    remove_outside_dfs(bond.neighbour, linker, keep_indices, visited, orig_visited);
  }
}

linker remove_outside_indices(linker const& linker, cif_data const& cif, double tolerance, std::string const& strategy)
{
  std::size_t starting_pos = 0;
  for(std::size_t i=0; i<linker.indices.size(); ++i)
  {
    if(linker.indices[i] < cif.num_atoms) { starting_pos = i; break; }
  }

  std::unordered_set<std::size_t> visited;
  std::unordered_set<std::size_t> orig_visited;
  std::vector<std::size_t> keep_indices;
  remove_outside_dfs(starting_pos, linker, keep_indices, visited, orig_visited);

  struct linker new_linker;
  for(std::size_t index : keep_indices)
  {
    new_linker.indices.push_back(linker.indices[index]);
  }
  new_linker.atoms = indices_to_atoms(new_linker.indices, cif.supercell);
  new_linker.graph = calc_connect_graph(new_linker.atoms, new_linker.atoms.size(), tolerance, strategy);
  new_linker.orig_indices = transform_to_original_indices(new_linker.indices, cif.num_atoms);
  new_linker.infinite = true;

  return new_linker;
}


std::string calc_cif_hash(cif_data& cif,
                          bool use_edge_attr, bool use_node_attr,
                          double tolerance,
                          std::string const& strategy)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1);
  cif.graph = calc_connect_graph(cif.supercell, cif.num_atoms, tolerance, strategy);
  label_graph(cif.graph, cif.num_atoms);
  auto reduced_graph = transform_to_original_indices(cif.graph, cif.num_atoms);
  auto hash = wl_graph_hash(reduced_graph, use_edge_attr, use_node_attr, 6, 16);

  return hash;
}

cif_data make_infinite_cif(linker const& linker, cif_data const& cif)
{
  auto new_cif = cif_data();
  new_cif.cell = cif.cell;
  new_cif.abs_atoms = linker.atoms;
  new_cif.num_atoms = linker.atoms.size();

  std::ofstream out("tmp.cif");
  print_cif(new_cif, &out);
  out.close();
  new_cif = cif_data("tmp.cif");
  return new_cif;
}

void decompose_infinite(cif_data& cif,
                      double tolerance,
                      std::string const& strategy, std::string const& base_filename)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1, false);
  cif.graph = calc_connect_graph_v2(cif.supercell, cif.supercell.size(), tolerance, strategy);


  std::vector<linker> linkers;
  std::vector<linker> nodes;

  std::unordered_set<std::size_t> visited;

  auto all_metal_indices = find_metal_idxs(cif.supercell);
  auto orig_cell_metal_indices = find_metal_idxs(cif.abs_atoms);
  if(all_metal_indices.empty()) { throw std::runtime_error("No metals found! Can't separate linker!"); }

  for(std::size_t metal_idx : all_metal_indices)
  {
    if(is_me_pillar(metal_idx, cif.graph, cif.supercell))
    {
      linker pillar;
      pillar.indices.push_back(metal_idx);
      for(auto const& bond : cif.graph.adjacency_list[metal_idx])
      {
        pillar.indices.push_back(bond.neighbour);
        visited.insert(bond.neighbour);
      }
      linkers.push_back(pillar);
    }
  }

	for(std::size_t metal_idx : all_metal_indices)
	{
    auto linker_starts = cif.graph.adjacency_list[metal_idx];
    for(auto const& bond : linker_starts)
    {
      if(visited.contains(bond.neighbour)) { continue; }
      if(data::is_metal(cif.supercell[bond.neighbour])) { continue; }
      std::vector<std::size_t> linker_indices;
      make_linker_dfs_v2(bond.neighbour, cif.graph, cif.supercell, linker_indices, visited);
      linkers.push_back({.indices=linker_indices});
    }
	}

  std::unordered_set<std::size_t> atoms_in_linkers;
  for(auto linker : linkers)
  {
    atoms_in_linkers.insert(std::begin(linker.indices), std::end(linker.indices));
  }

  std::unordered_set<std::size_t> visited2;
  for(std::size_t metal_idx : all_metal_indices)
  {
    std::vector<std::size_t> node_indices;
    if(visited2.contains(metal_idx)) { continue; }
    make_node_dfs_v2(metal_idx, cif.graph, cif.supercell, node_indices, visited2, atoms_in_linkers);
    nodes.push_back({.indices=node_indices});
  }

  std::unordered_set<std::size_t> atoms_in_nodes;
  for(auto node : nodes)
  {
    atoms_in_nodes.insert(std::begin(node.indices), std::end(node.indices));
  }


  remove_all_outside(linkers, cif.num_atoms);
  remove_all_outside(nodes, cif.num_atoms);

  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.graph = calc_connect_graph(linker.atoms, linker.atoms.size(), tolerance, strategy);
    linker.orig_indices = transform_to_original_indices(linker.indices, cif.num_atoms);
  }
  for(auto& node : nodes)
  {
    node.atoms = indices_to_atoms(node.indices, cif.supercell);
    node.graph = calc_connect_graph(node.atoms, node.atoms.size(), tolerance, strategy);
    node.orig_indices = transform_to_original_indices(node.indices, cif.num_atoms);
  }

  for(auto& linker : linkers)
  {
    linker.infinite = check_infinite(linker.indices, cif.num_atoms);
    if(linker.infinite)
    {
      linker = remove_outside_indices(linker, cif, tolerance, strategy);
      std::cerr << "WARNING! Infinite linker!\n";
    }
  }
  for(auto& node : nodes)
  {
    node.infinite = check_infinite(node.indices, cif.num_atoms);
    if(node.infinite)
    {
      node = remove_outside_indices(node, cif, tolerance, strategy);
      std::cerr << "WARNING! Infinite node!\n";
    }
  }


  remove_subsets(linkers);
  remove_subsets(nodes);

  // Get only numbers of unique linkers by HASH
  std::unordered_map<std::string, int> unique_linker_hashes;
  std::vector<std::size_t> unique_linker_indices;
  for(std::size_t i=0; i<linkers.size(); ++i)
  {
    linkers[i].graph.hash = wl_graph_hash(linkers[i].graph, false, true);
    if(unique_linker_hashes.contains(linkers[i].graph.hash))
    {
      unique_linker_hashes[linkers[i].graph.hash] += 1;
      continue;
    }
    unique_linker_hashes[linkers[i].graph.hash] = 1;
    unique_linker_indices.push_back(i);
  }

  std::size_t const n_unique_linkers = unique_linker_indices.size();
  std::vector<linker> unique_linkers(n_unique_linkers);
  for(std::size_t i=0; i<n_unique_linkers; ++i)
  {
    unique_linkers[i] = linkers[unique_linker_indices[i]];
  }


  // Get only numbers of unique nodes by HASH
  std::unordered_map<std::string, int> unique_node_hashes;
  std::vector<std::size_t> unique_node_indices;
  for(std::size_t i=0; i<nodes.size(); ++i)
  {
    nodes[i].graph.hash = wl_graph_hash(nodes[i].graph, false, true);
    if(unique_node_hashes.contains(nodes[i].graph.hash))
    {
      unique_node_hashes[nodes[i].graph.hash] += 1;
      continue;
    }
    unique_node_hashes[nodes[i].graph.hash] = 1;
    unique_node_indices.push_back(i);
  }

  std::size_t const n_unique_nodes = unique_node_indices.size();
  std::vector<linker> unique_nodes(n_unique_nodes);
  for(std::size_t i=0; i<n_unique_nodes; ++i)
  {
    unique_nodes[i] = nodes[unique_node_indices[i]];
  }

  // Put all neighbours in comments
  for(auto& linker : unique_linkers)
  {
    for(std::size_t i=0; i<linker.atoms.size(); ++i)
    {
      auto bonds = cif.graph.adjacency_list[linker.indices[i]];
      std::vector<atom> node_neighbours;
      for(auto const& bond : bonds)
      {
        if(atoms_in_nodes.contains(bond.neighbour)) { node_neighbours.push_back(cif.supercell[bond.neighbour]); }
      }
      if(node_neighbours.empty()) { continue; }
      linker.comment += std::to_string(i);
      linker.comment += ": ";
      linker.comment += std::to_string(node_neighbours.size());
      linker.comment += ", ";
      for(auto const& neighbour : node_neighbours)
      {
        linker.comment += std::to_string(neighbour.coord[0]);
        linker.comment += " ";
        linker.comment += std::to_string(neighbour.coord[1]);
        linker.comment += " ";
        linker.comment += std::to_string(neighbour.coord[2]);
        linker.comment += " ";
      }
      linker.comment += "; ";
    }
  }
  for(auto& node : unique_nodes)
  {
    for(std::size_t i=0; i<node.atoms.size(); ++i)
    {
      auto bonds = cif.graph.adjacency_list[node.indices[i]];
      std::vector<atom> linker_neighbours;
      for(auto const& bond : bonds)
      {
        if(atoms_in_linkers.contains(bond.neighbour)) { linker_neighbours.push_back(cif.supercell[bond.neighbour]); }
      }
      if(linker_neighbours.empty()) { continue; }
      node.comment += std::to_string(i);
      node.comment += ": ";
      node.comment += std::to_string(linker_neighbours.size());
      node.comment += ", ";
      for(auto const& neighbour : linker_neighbours)
      {
        node.comment += std::to_string(neighbour.coord[0]);
        node.comment += " ";
        node.comment += std::to_string(neighbour.coord[1]);
        node.comment += " ";
        node.comment += std::to_string(neighbour.coord[2]);
        node.comment += " ";
      }
      node.comment += "; ";
    }
  }

  // Write a file with numbers of building blocks
  std::ofstream ofs(base_filename + "_num_building_blocks.txt");
  for(std::size_t i=0; i<unique_linkers.size(); ++i)
  {
    ofs << "linker" << i;
    ofs << ".xyz: " << unique_linker_hashes[unique_linkers[i].graph.hash] << '\n';
  }
  for(std::size_t i=0; i<unique_nodes.size(); ++i)
  {
    ofs << "node" << i;
    ofs << ".xyz: " << unique_node_hashes[unique_nodes[i].graph.hash] << '\n';
  }

  for(std::size_t i=0; i<unique_linkers.size(); ++i)
  {
    std::string filename = base_filename;
    filename += "_linker";
    filename += std::to_string(i);
    filename += ".xyz";

    std::ofstream ofs(filename);
    print_xyz(unique_linkers[i].atoms, unique_linkers[i].comment, &ofs);
    print_graph_sane(unique_linkers[i].graph, &ofs);
  }
  for(std::size_t i=0; i<unique_nodes.size(); ++i)
  {
    std::string filename = base_filename;
    filename += "_node";
    filename += std::to_string(i);
    filename += ".xyz";

    std::ofstream ofs(filename);
    print_xyz(unique_nodes[i].atoms, unique_nodes[i].comment, &ofs);
    print_graph_sane(unique_nodes[i].graph, &ofs);
  }
}

void decompose_mof(cif_data& cif,
                   double tolerance,
                   std::string const& strategy,
                   bool use_linear_graph)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1);
  if(use_linear_graph)
  {
    cif.graph = calc_connect_graph_v3(cif.supercell, cif.fract_supercell, cif.supercell.size(), cif.cell, tolerance, strategy);
  }
  else
  {
    cif.graph = calc_connect_graph_v2(cif.supercell, cif.supercell.size(), tolerance, strategy);
  }

  std::vector<linker> linkers;
  std::vector<linker> nodes;
  std::vector<linker> disconnected;

  std::unordered_set<std::size_t> visited;

  auto all_metal_indices = find_metal_idxs(cif.supercell);
  auto orig_cell_metal_indices = find_metal_idxs(cif.abs_atoms);
  if(all_metal_indices.empty()) { throw std::runtime_error("No metals found! Can't separate linker!"); }

  for(std::size_t metal_idx : all_metal_indices)
  {
    if(is_me_pillar(metal_idx, cif.graph, cif.supercell))
    {
      linker pillar;
      pillar.indices.push_back(metal_idx);
      for(auto const& bond : cif.graph.adjacency_list[metal_idx])
      {
        pillar.indices.push_back(bond.neighbour);
        visited.insert(bond.neighbour);
      }
      linkers.push_back(pillar);
    }
  }

	for(std::size_t metal_idx : all_metal_indices)
	{
    auto linker_starts = cif.graph.adjacency_list[metal_idx];
    for(auto const& bond : linker_starts)
    {
      if(visited.contains(bond.neighbour)) { continue; }
      if(data::is_metal(cif.supercell[bond.neighbour])) { continue; }
      std::vector<std::size_t> linker_indices;
      make_linker_dfs_v2(bond.neighbour, cif.graph, cif.supercell, linker_indices, visited);
      if(linker_indices.size() > 4) { linkers.push_back({.indices=linker_indices}); }
    }
	}

  std::unordered_set<std::size_t> atoms_in_linkers;
  for(auto linker : linkers)
  {
    atoms_in_linkers.insert(std::begin(linker.indices), std::end(linker.indices));
  }

  std::unordered_set<std::size_t> visited2;
  for(std::size_t metal_idx : all_metal_indices)
  {
    std::vector<std::size_t> node_indices;
    if(visited2.contains(metal_idx)) { continue; }
    make_node_dfs_v2(metal_idx, cif.graph, cif.supercell, node_indices, visited2, atoms_in_linkers);
    nodes.push_back({.indices=node_indices});
  }

  // Drop atoms that are translated copies of linker atoms — these can leak in
  // from the supercell border.
  auto is_linker_translated = [&cif, &atoms_in_linkers](std::size_t index)
  {
    if(data::is_metal(cif.supercell[index])) { return false; }
    std::size_t reduced_index = reduce_index(index, cif.num_atoms);
    if(atoms_in_linkers.contains(reduced_index)) { return true; }
    return false;
  };

  for(auto& node : nodes)
  {
    auto new_end = std::remove_if(std::begin(node.indices), std::end(node.indices), is_linker_translated);
    node.indices.erase(new_end, std::end(node.indices));
  }

  // Add all atoms in nodes and into set
  std::unordered_set<std::size_t> atoms_in_nodes;
  for(auto node : nodes)
  {
    atoms_in_nodes.insert(std::begin(node.indices), std::end(node.indices));
  }

  // Find all disconnected with DFS
  std::unordered_set<std::size_t> visited3;
  for(std::size_t i=0; i<cif.supercell.size(); ++i)
  {
    if(atoms_in_linkers.contains(i) or atoms_in_nodes.contains(i)) { continue; }
    if(visited3.contains(i)) { continue; }
    std::vector<std::size_t> disconnected_indices;
    make_disconnected_dfs_v2(i, cif.graph, cif.supercell, disconnected_indices, visited3, atoms_in_linkers, atoms_in_nodes);
    disconnected.push_back({.indices=disconnected_indices});
  }


  remove_all_outside(linkers, cif.num_atoms);
  remove_all_outside(nodes, cif.num_atoms);
  remove_all_outside(disconnected, cif.num_atoms);

  for(auto& linker : linkers)
  {
    linker.atoms = indices_to_atoms(linker.indices, cif.supercell);
    linker.graph = calc_connect_graph(linker.atoms, linker.atoms.size(), tolerance, strategy);
    linker.orig_indices = transform_to_original_indices(linker.indices, cif.num_atoms);
  }
  for(auto& node : nodes)
  {
    node.atoms = indices_to_atoms(node.indices, cif.supercell);
    node.graph = calc_connect_graph(node.atoms, node.atoms.size(), tolerance, strategy);
    node.orig_indices = transform_to_original_indices(node.indices, cif.num_atoms);
  }
  for(auto& discon : disconnected)
  {
    discon.atoms = indices_to_atoms(discon.indices, cif.supercell);
    discon.graph = calc_connect_graph(discon.atoms, discon.atoms.size(), tolerance, strategy);
    discon.orig_indices = transform_to_original_indices(discon.indices, cif.num_atoms);
  }

  remove_subsets(linkers);
  remove_subsets(nodes);
  remove_subsets(disconnected);

  // Assign infinite status
  for(auto& linker : linkers)
  {
    linker.infinite = check_infinite(linker.indices, cif.num_atoms);
    if(linker.infinite)
    {
      linker = remove_outside_indices(linker, cif, tolerance, strategy);
      std::cerr << "WARNING! Infinite linker!\n";
    }
  }
  for(auto& node : nodes)
  {
    node.infinite = check_infinite(node.indices, cif.num_atoms);
    if(node.infinite)
    {
      node = remove_outside_indices(node, cif, tolerance, strategy);
      std::cerr << "WARNING! Infinite node!\n";
    }
  }
  for(auto& discon : disconnected)
  {
    discon.infinite = check_infinite(discon.indices, cif.num_atoms);
    if(discon.infinite)
    {
      discon = remove_outside_indices(discon, cif, tolerance, strategy);
      std::cerr << "WARNING! Infinite discon!\n";
    }
  }

  // Get only numbers of unique linkers by HASH
  std::unordered_map<std::string, std::vector<int>> unique_linker_hashes;
  std::vector<std::size_t> unique_linker_indices;
  for(std::size_t i=0; i<linkers.size(); ++i)
  {
    linkers[i].graph.hash = wl_graph_hash(linkers[i].graph, false, true);
    if(unique_linker_hashes.contains(linkers[i].graph.hash))
    {
      unique_linker_hashes[linkers[i].graph.hash].push_back(i);
      continue;
    }
    unique_linker_hashes[linkers[i].graph.hash] = std::vector<int> { static_cast<int>(i) };
    unique_linker_indices.push_back(i);
  }

  std::size_t const n_unique_linkers = unique_linker_indices.size();
  std::vector<linker> unique_linkers(n_unique_linkers);
  for(std::size_t i=0; i<n_unique_linkers; ++i)
  {
    unique_linkers[i] = linkers[unique_linker_indices[i]];
  }


  // Get only numbers of unique nodes by HASH
  std::unordered_map<std::string, std::vector<int>> unique_node_hashes;
  std::vector<std::size_t> unique_node_indices;
  for(std::size_t i=0; i<nodes.size(); ++i)
  {
    nodes[i].graph.hash = wl_graph_hash(nodes[i].graph, false, true);
    if(unique_node_hashes.contains(nodes[i].graph.hash))
    {
      unique_node_hashes[nodes[i].graph.hash].push_back(i);
      continue;
    }
    unique_node_hashes[nodes[i].graph.hash] = std::vector<int> { static_cast<int>(i) };
    unique_node_indices.push_back(i);
  }

  std::size_t const n_unique_nodes = unique_node_indices.size();
  std::vector<linker> unique_nodes(n_unique_nodes);
  for(std::size_t i=0; i<n_unique_nodes; ++i)
  {
    unique_nodes[i] = nodes[unique_node_indices[i]];
  }

  // Get only numbers of unique discons by HASH
  std::unordered_map<std::string, std::vector<int>> unique_discon_hashes;
  std::vector<std::size_t> unique_discon_indices;
  for(std::size_t i=0; i<disconnected.size(); ++i)
  {
    disconnected[i].graph.hash = wl_graph_hash(disconnected[i].graph, false, true);
    if(unique_discon_hashes.contains(disconnected[i].graph.hash))
    {
      unique_discon_hashes[disconnected[i].graph.hash].push_back(i);
      continue;
    }
    unique_discon_hashes[disconnected[i].graph.hash] = std::vector<int> { static_cast<int>(i) };
    unique_discon_indices.push_back(i);
  }

  std::size_t const n_unique_disconnected = unique_discon_indices.size();
  std::vector<linker> unique_disconnected(n_unique_disconnected);
  for(std::size_t i=0; i<n_unique_disconnected; ++i)
  {
    unique_disconnected[i] = disconnected[unique_discon_indices[i]];
  }

  // Put all neighbours in comments
  for(auto& linker : unique_linkers)
  {
    for(std::size_t i=0; i<linker.atoms.size(); ++i)
    {
      auto bonds = cif.graph.adjacency_list[linker.indices[i]];
      std::vector<atom> node_neighbours;
      for(auto const& bond : bonds)
      {
        if(atoms_in_nodes.contains(bond.neighbour)) { node_neighbours.push_back(cif.supercell[bond.neighbour]); }
      }
      if(node_neighbours.empty()) { continue; }
      linker.comment += std::to_string(i);
      linker.comment += ": ";
      linker.comment += std::to_string(node_neighbours.size());
      linker.comment += ", ";
      for(auto const& neighbour : node_neighbours)
      {
        linker.comment += std::to_string(neighbour.coord[0]);
        linker.comment += " ";
        linker.comment += std::to_string(neighbour.coord[1]);
        linker.comment += " ";
        linker.comment += std::to_string(neighbour.coord[2]);
        linker.comment += " ";
      }
      linker.comment += "; ";
    }
  }
  for(auto& node : unique_nodes)
  {
    for(std::size_t i=0; i<node.atoms.size(); ++i)
    {
      auto bonds = cif.graph.adjacency_list[node.indices[i]];
      std::vector<atom> linker_neighbours;
      for(auto const& bond : bonds)
      {
        if(atoms_in_linkers.contains(bond.neighbour)) { linker_neighbours.push_back(cif.supercell[bond.neighbour]); }
      }
      if(linker_neighbours.empty()) { continue; }
      node.comment += std::to_string(i);
      node.comment += ": ";
      node.comment += std::to_string(linker_neighbours.size());
      node.comment += ", ";
      for(auto const& neighbour : linker_neighbours)
      {
        node.comment += std::to_string(neighbour.coord[0]);
        node.comment += " ";
        node.comment += std::to_string(neighbour.coord[1]);
        node.comment += " ";
        node.comment += std::to_string(neighbour.coord[2]);
        node.comment += " ";
      }
      node.comment += "; ";
    }
  }

  // Write a file with numbers of building blocks
  std::ofstream ofs("num_building_blocks.txt");
  for(std::size_t i=0; i<unique_linkers.size(); ++i)
  {
    ofs << "linker" << i;
    if(unique_linkers[i].infinite) { ofs << "_infinite"; }
    ofs << ".xyz: " << unique_linker_hashes[unique_linkers[i].graph.hash].size() << '\n';
  }
  for(std::size_t i=0; i<unique_nodes.size(); ++i)
  {
    ofs << "node" << i;
    if(unique_nodes[i].infinite) { ofs << "_infinite"; }
    ofs << ".xyz: " << unique_node_hashes[unique_nodes[i].graph.hash].size() << '\n';
  }
  for(std::size_t i=0; i<unique_disconnected.size(); ++i)
  {
    ofs << "disconnected" << i;
    if(unique_disconnected[i].infinite) { ofs << "_infinite"; }
    ofs << ".xyz: " << unique_discon_hashes[unique_disconnected[i].graph.hash].size() << '\n';
  }

  for(std::size_t i=0; i<unique_linkers.size(); ++i)
  {
    std::string filename = "linker";
    filename += std::to_string(i);
    if(unique_linkers[i].infinite) { filename += "_infinite"; }
    filename += ".xyz";

    std::ofstream ofs(filename);
    print_xyz(unique_linkers[i].atoms, unique_linkers[i].comment, &ofs);
    print_graph_sane(unique_linkers[i].graph, &ofs);
  }
  for(std::size_t i=0; i<unique_nodes.size(); ++i)
  {
    std::string filename = "node";
    filename += std::to_string(i);
    if(unique_nodes[i].infinite)
    {
      filename += "_infinite";
      auto new_cif = make_infinite_cif(unique_nodes[i], cif);
      decompose_infinite(new_cif, tolerance, strategy, filename);
    }
    filename += ".xyz";

    std::ofstream ofs(filename);
    print_xyz(unique_nodes[i].atoms, unique_nodes[i].comment, &ofs);
    print_graph_sane(unique_nodes[i].graph, &ofs);
  }
  for(std::size_t i=0; i<unique_disconnected.size(); ++i)
  {
    std::string filename = "disconnected";
    filename += std::to_string(i);
    if(unique_disconnected[i].infinite) { filename += "_infinite"; }
    filename += ".xyz";

    std::ofstream ofs2(filename);
    print_xyz(unique_disconnected[i].atoms, "", &ofs2);
    print_graph_sane(unique_disconnected[i].graph, &ofs2);
  }


  // DUMP all non-unique structures for everything
  for(std::size_t i=0; i<unique_linkers.size(); ++i)
  {
    std::string filename = "DUMP_linker";
    filename += std::to_string(i);
    if(unique_linkers[i].infinite) { filename += "_infinite"; }
    filename += ".xyz";

    std::ofstream ofs(filename);
    for(auto const& non_unique_linkers_index : unique_linker_hashes.at(unique_linkers[i].graph.hash))
    {
      print_xyz(linkers[non_unique_linkers_index].atoms, "NON-UNIQUE-CONFORMATION", &ofs);
      ofs << '\n';
    }
  }
  for(std::size_t i=0; i<unique_nodes.size(); ++i)
  {
    std::string filename = "DUMP_node";
    filename += std::to_string(i);
    if(unique_nodes[i].infinite) { filename += "_infinite"; }
    filename += ".xyz";

    std::ofstream ofs(filename);
    for(auto const& non_unique_nodes_index : unique_node_hashes.at(unique_nodes[i].graph.hash))
    {
      print_xyz(nodes[non_unique_nodes_index].atoms, "NON-UNIQUE-CONFORMATION", &ofs);
      ofs << '\n';
    }
  }
  for(std::size_t i=0; i<unique_disconnected.size(); ++i)
  {
    std::string filename = "DUMP_disconnected";
    filename += std::to_string(i);
    if(unique_disconnected[i].infinite) { filename += "_infinite"; }
    filename += ".xyz";

    std::ofstream ofs(filename);
    for(auto const& non_unique_discon_index : unique_discon_hashes.at(unique_disconnected[i].graph.hash))
    {
      print_xyz(disconnected[non_unique_discon_index].atoms, "NON-UNIQUE-CONFORMATION", &ofs);
      ofs << '\n';
    }
  }
}

#endif
