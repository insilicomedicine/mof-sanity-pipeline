// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "graph_ops.hpp"
#include "cif.hpp"
#include "cif_ops.hpp"
#include "data_tables.hpp"
#include "linalg.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <sstream>

std::vector<bond> find_bonds(std::size_t target, std::vector<atom> const& atoms,
                             double tolerance, std::string const& strategy)
{
  std::vector<bond> neighbours;
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    if(target==i) { continue; }
    double const r = atoms_range(atoms[target], atoms[i]);
    if(data::is_bonded(atoms[target].name, atoms[i].name, r, strategy, tolerance))
    {
      if(data::is_metal(atoms[target]) or data::is_metal(atoms[i]))
      {
        neighbours.push_back({i, r, bond_type::METAL});
      }
      else
      {
        neighbours.push_back({i,r, bond_type::SINGLE});
      }
    }
  }
  return neighbours;
}

molecular_graph calc_connect_graph(std::vector<atom> const& atoms, std::size_t graph_size,
                                   double tolerance, std::string const& strategy)
{
  molecular_graph graph;
  graph.adjacency_list.resize(graph_size);
  graph.node_labels.resize(graph_size);

	for(std::size_t i=0; i<graph_size; ++i)
	{
    graph.adjacency_list[i] = find_bonds(i, atoms, tolerance, strategy);
    graph.node_labels[i] = atoms[i].name;
	}

  return graph;
}

molecular_graph calc_connect_graph_v2(std::vector<atom> const& atoms, std::size_t graph_size,
                                      double tolerance, std::string const& strategy)
{
  auto ranges_matrix = std::vector<std::vector<float>>(atoms.size(), std::vector<float>(atoms.size()));

  // Upper left block
  for(std::size_t i=0; i<graph_size; ++i)
  {
    for(std::size_t j=i+1; j<graph_size; ++j)
    {
      double r = atoms_range(atoms[i], atoms[j]);
      ranges_matrix[i][j] = r;
      ranges_matrix[j][i] = r;
    }
  }

  // Upper right block
  for(std::size_t i=0; i<graph_size; ++i)
  {
    for(std::size_t j=graph_size; j<atoms.size(); ++j)
    {
      double r = atoms_range(atoms[i], atoms[j]);
      ranges_matrix[i][j] = r;
    }
  }

  molecular_graph graph;
  graph.adjacency_list.resize(graph_size);
  graph.node_labels.resize(graph_size);

  for(std::size_t i=0; i<graph_size; ++i) { graph.node_labels[i] = atoms[i].name; }

  std::vector<double> cov_radii(graph_size);
  for(std::size_t i=0; i<graph_size; ++i)
  {
    cov_radii[i] = data::get_cov_radius(atoms[i].name, strategy);
  }


  for(std::size_t i=0; i<graph_size; ++i)
  {
    for(std::size_t j=i+1; j<graph_size; ++j)
    {
      double r = ranges_matrix[i][j];
      double cov_rad_i = cov_radii[i];
      double cov_rad_j = cov_radii[j];

      if(strategy=="jmol")
      {
        if((atoms[i].name=="H" or atoms[j].name=="H") and
           (atoms[i].name=="P" or atoms[j].name=="P")) { cov_rad_i += 1.0; }
      }

      if(data::is_bonded(cov_rad_i, cov_rad_j, r, tolerance))
      {
        auto bt = (data::is_metal(atoms[i]) or data::is_metal(atoms[j])) ? bond_type::METAL : bond_type::SINGLE;
        graph.adjacency_list[i].push_back({j, r, bt});
        graph.adjacency_list[j].push_back({i, r, bt});
      }
    }
  }

  for(std::size_t i=0; i<graph_size; ++i)
  {
    for(std::size_t j=graph_size; j<atoms.size(); ++j)
    {
      double r = ranges_matrix[i][j];
      if(data::is_bonded(atoms[i].name, atoms[j].name, r, strategy, tolerance))
      {
        auto bt = (data::is_metal(atoms[i]) or data::is_metal(atoms[j])) ? bond_type::METAL : bond_type::SINGLE;
        graph.adjacency_list[i].push_back({j, r, bt});
      }
    }
  }

  return graph;
}

double max_bond_length(std::vector<double> const& cov_radii, double tolerance)
{
  double max_radius = *std::max_element(std::begin(cov_radii), std::end(cov_radii));
  return data::apply_tolerance(2.0*max_radius, tolerance) + 0.05;
}

struct grid
{
  int size_a;
  int size_b;
  int size_c;
  double len_a;
  double len_b;
  double len_c;
  std::vector<std::vector<std::vector<std::vector<std::size_t>>>> cells;

  grid(std::vector<atom> const& atoms, unit_cell const& cell, double max_bond_length, double skin = 0.4)
  {
    double len_a_supercell = 3.0*norm(vec3d(cell.axes.a));
    double len_b_supercell = 3.0*norm(vec3d(cell.axes.b));
    double len_c_supercell = 3.0*norm(vec3d(cell.axes.c));

    size_a = static_cast<int>(std::ceil(len_a_supercell / (max_bond_length+skin)));
    size_b = static_cast<int>(std::ceil(len_b_supercell / (max_bond_length+skin)));
    size_c = static_cast<int>(std::ceil(len_c_supercell / (max_bond_length+skin)));
    len_a = 3.0 / size_a;
    len_b = 3.0 / size_b;
    len_c = 3.0 / size_c;

    cells.resize(size_a);
    for(auto& item : cells) { item.resize(size_b); }
    for(auto& i : cells)
    {
      for(auto& j : i) { j.resize(size_c); }
    }

    for(std::size_t i=0; i<atoms.size(); ++i)
    {
      int idx_a = static_cast<int>(std::floor((atoms[i].coord[0]+1) / len_a));
      int idx_b = static_cast<int>(std::floor((atoms[i].coord[1]+1) / len_b));
      int idx_c = static_cast<int>(std::floor((atoms[i].coord[2]+1) / len_c));
      if(idx_a<0 or idx_b<0 or idx_c<0) { print_atom(atoms[i], &std::cerr); throw std::runtime_error("INDEX LESS THAN 0!"); }
      if(idx_a>=size_a or idx_b>=size_b or idx_c>=size_c) { print_atom(atoms[i], &std::cerr); throw std::runtime_error("INDEX BIGGER THAN SIZE!"); }
      cells[idx_a][idx_b][idx_c].push_back(i);
    }

  }
};

molecular_graph calc_connect_graph_v3(std::vector<atom> const& atoms, std::vector<atom> const& fract_atoms,
                                      std::size_t graph_size, unit_cell const& cell,
                                      double tolerance, std::string const& strategy)
{
  molecular_graph graph;
  graph.adjacency_list.resize(graph_size);
  graph.node_labels.resize(graph_size);

  for(std::size_t i=0; i<graph_size; ++i) { graph.node_labels[i] = atoms[i].name; }

  std::vector<double> cov_radii(atoms.size());
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    cov_radii[i] = data::get_cov_radius(atoms[i].name, strategy);
  }

  double max_bond = max_bond_length(cov_radii, tolerance);
  grid g(fract_atoms, cell, max_bond);

  constexpr double skin = 0.4;
  auto const compute_max_offset = [&](double len_supercell, int num_cells)
  {
    if(num_cells <= 0) { return 0; }
    double const cell_length = len_supercell / static_cast<double>(num_cells);
    if(cell_length <= 0.0) { return num_cells; }
    double const ratio = (max_bond + skin) / cell_length;
    int const offset = static_cast<int>(std::ceil(ratio));
    return std::max(1, offset);
  };

  double const len_a_supercell = 3.0 * norm(vec3d(cell.axes.a));
  double const len_b_supercell = 3.0 * norm(vec3d(cell.axes.b));
  double const len_c_supercell = 3.0 * norm(vec3d(cell.axes.c));

  int const max_offset_a = compute_max_offset(len_a_supercell, g.size_a);
  int const max_offset_b = compute_max_offset(len_b_supercell, g.size_b);
  int const max_offset_c = compute_max_offset(len_c_supercell, g.size_c);

  // Match calc_connect_graph_v2 semantics: store bonds only once and limit
  // adjacency entries to atoms inside the requested subgraph.
  auto try_add_edge = [&](std::size_t idx1, std::size_t idx2)
  {
    if(idx2 >= idx1) { return; }

    bool const idx1_in_graph = idx1 < graph_size;
    bool const idx2_in_graph = idx2 < graph_size;
    if(!idx1_in_graph && !idx2_in_graph) { return; }

    double const r = atoms_range(atoms[idx1], atoms[idx2]);

    double cov_rad1 = cov_radii[idx1];
    double cov_rad2 = cov_radii[idx2];
    if(strategy=="jmol")
    {
      if((atoms[idx1].name=="H" or atoms[idx2].name=="H") and
          (atoms[idx1].name=="P" or atoms[idx2].name=="P")) { cov_rad1 += 1.0; }
    }
    if(!data::is_bonded(cov_rad1, cov_rad2, r, tolerance)) { return; }

    auto const bt = (data::is_metal(atoms[idx1]) or data::is_metal(atoms[idx2])) ? bond_type::METAL : bond_type::SINGLE;

    if(idx1_in_graph)
    {
      graph.adjacency_list[idx1].push_back({idx2, r, bt});
    }
    if(idx2_in_graph)
    {
      graph.adjacency_list[idx2].push_back({idx1, r, bt});
    }
  };

  for(int i=0; i<g.size_a; ++i)
  {
    for(int j=0; j<g.size_b; ++j)
    {
      for(int k=0; k<g.size_c; ++k)
      {
        // Inside cell
        for(std::size_t idx1 : g.cells[i][j][k])
        {
          for(std::size_t idx2 : g.cells[i][j][k])
          {
            try_add_edge(idx1, idx2);
          }
        }
        // Neighbouring cells
        for(int di=-max_offset_a; di<=max_offset_a; ++di)
        {
          int new_i = i + di;
          if(new_i<0 or new_i >= g.size_a) { continue; }
          for(int dj=-max_offset_b; dj<=max_offset_b; ++dj)
          {
            int new_j = j + dj;
            if(new_j<0 or new_j >= g.size_b) { continue; }
            for(int dk=-max_offset_c; dk<=max_offset_c; ++dk)
            {
              int new_k = k + dk;
              if(new_k<0 or new_k>=g.size_c) { continue; }
              if(di==0 and dj==0 and dk==0) { continue; }

              for(std::size_t idx1 : g.cells[i][j][k])
              {
                for(std::size_t idx2 : g.cells[new_i][new_j][new_k])
                {
                  try_add_edge(idx1, idx2);
                }
              }
            }
          }
        }
      }
    }
  }

  for(auto& neighbours : graph.adjacency_list)
  {
    std::sort(std::begin(neighbours), std::end(neighbours),
              [](bond const& lhs, bond const& rhs)
              {
                if(lhs.neighbour != rhs.neighbour) { return lhs.neighbour < rhs.neighbour; }
                return lhs.r < rhs.r;
              });
  }
  return graph;
}


bool check_overlap(molecular_graph const& graph, std::vector<atom> const& atoms, double r_min)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    for(auto const& bond : graph.adjacency_list[i])
    {
      if(bond.r < r_min)
      {
        std::cerr << "ERROR!! Atoms " << atoms[i].name << " #" << i+1
                  << " and " << atoms[bond.neighbour].name << " #" << bond.neighbour+1
                  << " are closer than " << r_min << std::endl;
        return false;
      }
    }
  }
  return true;
}

void remove_too_long_bonds(molecular_graph& graph, double r_max)
{
  for(auto& bonds : graph.adjacency_list)
  {
      bonds.erase(std::remove_if(std::begin(bonds), std::end(bonds),
                                 [r_max](auto const& b){ return b.r > r_max; }),
                  std::end(bonds));
  }
}

bool calc_graph(cif_data& cif,
                double r_min, double r_max,
                double tolerance,
                std::string const& strategy)
{
  cif.calc_supercell(-1, 1, -1, 1, -1, 1);
  cif.graph = calc_connect_graph(cif.supercell, cif.num_atoms, tolerance, strategy);
  if(!check_overlap(cif.graph, cif.supercell, r_min)) { return false; }
  remove_too_long_bonds(cif.graph, r_max);
  return true;
}

bool verify_trim_graph(molecular_graph &graph, std::vector<atom> const& atoms, double r_min, double r_max)
{
  if(!check_overlap(graph, atoms, r_min)) { return false; }
  remove_too_long_bonds(graph, r_max);
  return true;
}

bool all_atoms_unsaturated(std::vector<std::size_t> const& targets,
                           std::vector<atom> const& atoms, molecular_graph const& graph)
{
  for(std::size_t target : targets)
  {
    std::size_t curr_valence = graph.adjacency_list[target].size();
    if(!data::is_unsaturated(atoms[target].name, curr_valence)) { return false; }
  }
  return true;
}

void find_cycles_dfs(std::size_t n, std::size_t target,
                     std::vector<std::vector<std::size_t>>& rings, molecular_graph const& graph,
                     std::vector<int>& colors, std::vector<std::size_t>& traversal_order, int parent = -1)
{
  colors[target] = 1;
  traversal_order.push_back(target);

  for(auto const& bond : graph.adjacency_list[target])
  {
    if(colors[bond.neighbour] == 0)
    {
      if(graph.adjacency_list[bond.neighbour].size() < 2) { colors[bond.neighbour] = 2; }
      else
      {
        find_cycles_dfs(n, bond.neighbour, rings, graph, colors, traversal_order, target);
      }
    }
    else if(colors[bond.neighbour] == 1)
    {
      if(parent != -1 and bond.neighbour != static_cast<std::size_t>(parent))
      {
        std::vector<std::size_t> ring;
        auto last_it = std::find(traversal_order.rbegin(), traversal_order.rend(), target);
        for(auto rit = last_it; rit!=traversal_order.rend() && (*rit)!=bond.neighbour; ++rit)
        {
          ring.push_back(*rit);
        }
        ring.push_back(bond.neighbour);
        if(ring.size() == 6) { rings.push_back(ring); }
      }
    }
  }

  colors[target] = 2;
  traversal_order.pop_back();
}

std::vector<std::vector<std::size_t>> find_all_n_member_rings(molecular_graph const& graph, std::size_t n)
{
  std::vector<std::vector<size_t>> rings;
  std::size_t const num_atoms = graph.adjacency_list.size();
  std::vector<int> colors(num_atoms, 0);

  for(std::size_t i=0; i<num_atoms; ++i)
  {
    if(colors[i] != 0) { continue; }
    if(graph.adjacency_list[i].size() < 2) { colors[i] = 2; }
    std::vector<std::size_t> traversal_order;
    find_cycles_dfs(n, i, rings, graph, colors, traversal_order);
  }

  return rings;
}

void make_linker_dfs_v2(std::size_t target,
                        molecular_graph const& graph,
                        std::vector<atom> const& atoms,
                        std::vector<std::size_t>& linker,
                        std::unordered_set<std::size_t>& visited)
{
  if(data::is_metal(atoms[target])) { return; }
	linker.push_back(target);
	visited.insert(target);

  for(auto const& bond : graph.adjacency_list[target])
  {
    if(visited.contains(bond.neighbour)) { continue; }
    make_linker_dfs_v2(bond.neighbour, graph, atoms, linker, visited);
  }
}

void make_node_dfs_v2(std::size_t target,
                      molecular_graph const& graph,
                      std::vector<atom> const& atoms,
                      std::vector<std::size_t>& node,
                      std::unordered_set<std::size_t>& visited,
                      std::unordered_set<std::size_t> const& linkers_indices)
{
  if(linkers_indices.contains(target)) { return; }
  node.push_back(target);
  visited.insert(target);

  for(auto const& bond : graph.adjacency_list[target])
  {
    if(visited.contains(bond.neighbour)) { continue; }
    make_node_dfs_v2(bond.neighbour, graph, atoms, node, visited, linkers_indices);
  }
}

void make_disconnected_dfs_v2(std::size_t target,
                      molecular_graph const& graph,
                      std::vector<atom> const& atoms,
                      std::vector<std::size_t>& node,
                      std::unordered_set<std::size_t>& visited,
                      std::unordered_set<std::size_t> const& linkers_indices,
                      std::unordered_set<std::size_t> const& nodes_indices)
{
  if(linkers_indices.contains(target)) { return; }
  if(nodes_indices.contains(target)) { return; }
  node.push_back(target);
  visited.insert(target);

  for(auto const& bond : graph.adjacency_list[target])
  {
    if(visited.contains(bond.neighbour)) { continue; }
    make_disconnected_dfs_v2(bond.neighbour, graph, atoms, node, visited, linkers_indices, nodes_indices);
  }
}

std::vector<std::size_t> find_overlap_atoms(molecular_graph const& graph, double r_min)
{
  std::vector<std::size_t> overlap_indices;
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    for(auto const& bond : graph.adjacency_list[i])
    {
      if(bond.r < r_min)
      {
        overlap_indices.push_back(i);
      }
    }
  }
  return overlap_indices;
}

void find_connected_dfs(std::size_t target, std::unordered_set<std::size_t>& visited, std::vector<std::size_t>& subgraph, molecular_graph const& graph)
{
  visited.insert(target);
  subgraph.push_back(target);
  for(auto const& bond : graph.adjacency_list[target])
  {
    if(visited.contains(bond.neighbour)) { continue; }
    find_connected_dfs(bond.neighbour, visited, subgraph, graph);
  }
}

molecular_graph transform_to_original_indices(molecular_graph graph, std::size_t num_atoms)
{
  for(auto& bonds : graph.adjacency_list)
  {
    for(auto& bond : bonds) { bond.neighbour = bond.neighbour % num_atoms; }
  }
  return graph;
}

bool check_disconnected(cif_data const& cif, std::size_t max_err_size)
{
  auto reduced_graph = transform_to_original_indices(cif.graph, cif.num_atoms);
  std::unordered_set<std::size_t> visited;
  std::vector<std::vector<std::size_t>> subgraphs;
  for(std::size_t i=0; i<cif.num_atoms; ++i)
  {
    if(visited.contains(i)) { continue; }
    std::vector<std::size_t> subgraph;
    find_connected_dfs(i, visited, subgraph, reduced_graph);
    subgraphs.push_back(subgraph);
  }

  if(subgraphs.size() == 1) { return true; }

  // EDGE CASE #1: DO NOT count isolated subgraphs only consisting of metals
  auto no_metals_end = std::remove_if(std::begin(subgraphs), std::end(subgraphs),
                                      [&cif](auto const& vec)
                                      {
                                        return std::all_of(std::begin(vec), std::end(vec), [&cif](std::size_t n)
                                          {
                                            return data::is_metal(cif.supercell[n]);
                                          });
                                      });
  subgraphs.erase(no_metals_end, std::end(subgraphs));
  // If valid without metals
  if(subgraphs.size() == 1) { return true; }

  // EDGE CASE #2: isolated single atom (special message)
  for(auto const& subgraph : subgraphs)
  {
    if(subgraph.size() == 1)
    {
      std::cerr << "ERROR!! Atom " << cif.abs_atoms[subgraph[0]].name << " #" << subgraph[0]+1 << " has NO neighbours" << std::endl;
      return false;
    }
  }

  // EDGE CASE #3: isolated subgraph connected to other cell in supercell
  std::vector<std::unordered_set<std::size_t>> uniq_graphs;
  for(auto const& subgraph : subgraphs)
  {
    auto has_elem_set_it = std::end(uniq_graphs); // not end if has any same element

    for(std::size_t num : subgraph)
    {
      has_elem_set_it = std::find_if(std::begin(uniq_graphs), std::end(uniq_graphs),
                                      [num](auto&& uniq_graph){ return uniq_graph.contains(num); });
      if(has_elem_set_it != std::end(uniq_graphs)) { break; }
    }

    if(has_elem_set_it == std::end(uniq_graphs))
    {
      uniq_graphs.emplace_back(std::begin(subgraph), std::end(subgraph));
    }
    else
    {
      has_elem_set_it->insert(std::begin(subgraph), std::end(subgraph));
    }
  }

  if(uniq_graphs.size() > 1)
  {
	  for(auto const& uniq_graph : uniq_graphs)
	  {
		  if(uniq_graph.size() > max_err_size)
		  {
			  continue;
		  }
		  else
		  {
        std::cerr << "ERROR!! Graph contains disconnected subgraphs!" << std::endl;
        for(std::size_t n=0; n<subgraphs.size(); ++n)
        {
          std::cerr << "SUBGRAPH #" << n+1 << std::endl;
          for(std::size_t num : subgraphs[n]) { std::cerr << num + 1 << ' '; }
          std::cerr << std::endl;
          auto subgraph_final_atoms = indices_to_atoms(subgraphs[n], cif.supercell);
          print_xyz(subgraph_final_atoms, "", &std::cerr);
          std::cerr << std::endl;
          std::cerr << std::endl;
        }
        return false;
		  }
	  }
  }

  return true;
}

std::vector<std::vector<std::size_t>> find_subgraphs(molecular_graph const& graph)
{
  std::unordered_set<std::size_t> visited;
  std::vector<std::vector<std::size_t>> subgraphs;
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    if(visited.contains(i)) { continue; }
    else { std::cerr << "visited has no number " << i << std::endl; }
    std::vector<std::size_t> subgraph;
    find_connected_dfs(i, visited, subgraph, graph);
    subgraphs.push_back(subgraph);
  }
  return subgraphs;
}

std::string carbon_hybridization(std::size_t num_neighbours)
{
  std::string hybridization;
  switch(num_neighbours)
  {
    case 4:
      hybridization = "sp3"; break;
    case 3:
      hybridization = "sp2"; break;
    case 2:
      hybridization = "sp"; break;
    case 1:
      hybridization = "sp"; break;
  }
  return hybridization;
}

enum class N_HIBRIDIZATION
{
  SP3, SP2, SP, NONE
};

std::vector<N_HIBRIDIZATION> nitrogen_hybridization(std::size_t target, molecular_graph const& graph, [[maybe_unused]] std::vector<atom> const& atoms)
{
  std::size_t num_neighbours = graph.adjacency_list[target].size();
  switch(num_neighbours)
  {
    case 4:
      return {N_HIBRIDIZATION::SP3};
    case 3:
      return {N_HIBRIDIZATION::SP3, N_HIBRIDIZATION::SP2};
    case 2:
      return {N_HIBRIDIZATION::SP2, N_HIBRIDIZATION::SP};
    case 1:
      return {N_HIBRIDIZATION::SP};
    default:
      return {N_HIBRIDIZATION::NONE};
  }
}

bool check_carbon_sp(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  if(bonds.size() == 1) { return true; }
  if(bonds.size() != 2) { throw std::runtime_error("Wrong atom provided for sp check!"); }

  constexpr double sp_angle_proper = ang_to_rad(180);
  constexpr double sp_angle_tolerance = ang_to_rad(10);

  auto c = vec3d(atoms[target].coord);
  auto a = vec3d(atoms[bonds[0].neighbour].coord);
  auto b = vec3d(atoms[bonds[1].neighbour].coord);

  auto ca = a - c;
  auto cb = b - c;

  double theta = angle(ca, cb);

  if(std::abs(theta - sp_angle_proper) > sp_angle_tolerance)
  {
    std::cerr << "BAD sp angle detected for atom C #" << target+1 << std::endl;
    return false;
  }
  return true;
}

bool check_nitrogen_sp(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  if(bonds.size() == 1) { return true; }
  if(bonds.size() != 2) { throw std::runtime_error("Wrong atom provided for sp check!"); }

  constexpr double sp_angle_proper = ang_to_rad(180);
  constexpr double sp_angle_tolerance = ang_to_rad(10);

  auto c = vec3d(atoms[target].coord);
  auto a = vec3d(atoms[bonds[0].neighbour].coord);
  auto b = vec3d(atoms[bonds[1].neighbour].coord);

  auto ca = a - c;
  auto cb = b - c;

  double theta = angle(ca, cb);

  if(std::abs(theta - sp_angle_proper) > sp_angle_tolerance)
  {
    std::cerr << "BAD sp angle detected for atom N #" << target+1 << std::endl;
    std::cerr << "Proper = " << sp_angle_proper << " Actual = " << theta << std::endl;
    return false;
  }
  return true;
}

bool check_carbon_sp2(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  if(bonds.size() != 3) { throw std::runtime_error("Wrong atom provided for sp2 check!"); }

  constexpr double sp2_angle_proper = ang_to_rad(120);
  constexpr double sp2_angle_tolerance = ang_to_rad(10);
  constexpr double sp2_dihedral_proper0 = ang_to_rad(0);
  constexpr double sp2_dihedral_proper180 = ang_to_rad(180);
  constexpr double sp2_dihedral_tolerance = ang_to_rad(10);

  auto x = vec3d(atoms[target].coord);
  auto a = vec3d(atoms[bonds[0].neighbour].coord);
  auto b = vec3d(atoms[bonds[1].neighbour].coord);
  auto c = vec3d(atoms[bonds[2].neighbour].coord);

  auto xa = a - x;
  auto xb = b - x;
  auto xc = c - x;

  auto normal = vector_product(xa, xb);
  double theta = M_PI_2 - angle(normal, xc);

  if(!within_range(theta, sp2_dihedral_proper0, sp2_angle_tolerance) and
     !within_range(theta, sp2_dihedral_proper180, sp2_dihedral_tolerance))
  {
    std::cerr << "BAD sp2 dihedral detected for atom C #" << target+1 << std::endl;
    return false;
  }

  double axb = angle(xa, xb);
  double axc = angle(xa, xc);
  double bxc = angle(xb, xc);

  if(!within_range(axb, sp2_angle_proper, sp2_angle_tolerance) or
     !within_range(axc, sp2_angle_proper, sp2_angle_tolerance) or
     !within_range(bxc, sp2_angle_proper, sp2_angle_tolerance))
  {
    std::cerr << "BAD sp2 angle detected for atom C #" << target+1 << std::endl;
    return false;
  }

  return true;
}

bool check_nitrogen_sp2(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  constexpr double sp2_angle_proper = ang_to_rad(120);
  constexpr double sp2_angle_tolerance = ang_to_rad(10);

  auto x = vec3d(atoms[target].coord);
  std::vector<vec3d> neighbours;
  for(auto const& bond : bonds)
  {
    neighbours.emplace_back(atoms[bond.neighbour].coord);
  }

  std::vector<vec3d> ranges;
  for(auto const& neighbour : neighbours)
  {
    ranges.push_back(neighbour - x);
  }

  std::vector<double> angles;
  for(std::size_t i=0; i<ranges.size(); ++i)
  {
    for(std::size_t j=i+1; j<ranges.size(); ++j)
    {
      angles.push_back(angle(ranges[i], ranges[j]));
    }
  }

  for(double angle : angles)
  {
    if(!within_range(angle, sp2_angle_proper, sp2_angle_tolerance))
    {
      std::cerr << "BAD sp2 angle detected for atom N #" << target+1 << std::endl;
      std::cerr << "Proper = " << sp2_angle_proper << " Actual = " << angle << std::endl;
      return false;
    }
  }

  return true;
}

bool check_carbon_sp3(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  if(bonds.size() != 4) { throw std::runtime_error("Wrong atom provided for sp3 check!"); }

  constexpr double sp3_angle_proper = ang_to_rad(109);
  constexpr double sp3_angle_tolerance = ang_to_rad(10);

  auto x = vec3d(atoms[target].coord);
  auto a = vec3d(atoms[bonds[0].neighbour].coord);
  auto b = vec3d(atoms[bonds[1].neighbour].coord);
  auto c = vec3d(atoms[bonds[2].neighbour].coord);
  auto d = vec3d(atoms[bonds[3].neighbour].coord);

  auto xa = a - x;
  auto xb = b - x;
  auto xc = c - x;
  auto xd = d - x;
  
  double axb = angle(xa, xb);
  double axc = angle(xa, xc);
  double axd = angle(xa, xd);
  double bxc = angle(xb, xc);
  double bxd = angle(xb, xd);
  double cxd = angle(xc, xd);

  if(!within_range(axb, sp3_angle_proper, sp3_angle_tolerance) or
     !within_range(axc, sp3_angle_proper, sp3_angle_tolerance) or
     !within_range(axd, sp3_angle_proper, sp3_angle_tolerance) or
     !within_range(bxc, sp3_angle_proper, sp3_angle_tolerance) or
     !within_range(bxd, sp3_angle_proper, sp3_angle_tolerance) or
     !within_range(cxd, sp3_angle_proper, sp3_angle_tolerance)
  )
  {
    std::cerr << "BAD sp3 angle detected for atom C #" << target+1 << std::endl;
    return false;
  }
  return true;
}

bool check_nitrogen_sp3(std::size_t target, std::vector<bond> const& bonds, std::vector<atom> const& atoms)
{
  constexpr double sp3_angle_proper = ang_to_rad(107);
  constexpr double sp3_angle_tolerance = ang_to_rad(10);

  auto x = vec3d(atoms[target].coord);
  std::vector<vec3d> neighbours;
  for(auto const& bond : bonds)
  {
    neighbours.emplace_back(atoms[bond.neighbour].coord);
  }

  std::vector<vec3d> ranges(neighbours.size());
  for(std::size_t i=0; i<neighbours.size(); ++i)
  {
    ranges[i] = neighbours[i] - x;
  }

  std::vector<double> angles;
  for(std::size_t i=0; i<ranges.size(); ++i)
  {
    for(std::size_t j=i+1; j<ranges.size(); ++j)
    {
      angles.push_back(angle(ranges[i], ranges[j]));
    }
  }

  for(double angle : angles)
  {
    if(!within_range(angle, sp3_angle_proper, sp3_angle_tolerance))
    {
      std::cerr << "BAD sp3 angle detected for atom N #" << target+1 << std::endl;
      std::cerr << "Proper = " << sp3_angle_proper << " Actual = " << angle << std::endl;
      return false;
    }
  }
  return true;
}

bool check_carbons(molecular_graph const& graph, std::vector<atom> const& atoms)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    if(atoms[i].name != "C") { continue; }
    std::string hybridization = carbon_hybridization(graph.adjacency_list[i].size());
    if(hybridization == "sp")
    {
      if(!check_carbon_sp(i, graph.adjacency_list[i], atoms)) { return false; }
    }
    if(hybridization == "sp2")
    {
      if(!check_carbon_sp2(i, graph.adjacency_list[i], atoms)) { return false; }
    }
    if(hybridization == "sp3")
    {
      if(!check_carbon_sp3(i, graph.adjacency_list[i], atoms)) { return false; }
    }
  }
  return true;
}

bool check_nitrogens(molecular_graph const& graph, std::vector<atom> const& atoms)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    if(atoms[i].name != "N") { continue; }

    std::size_t num_neighbours = graph.adjacency_list[i].size();
    if(num_neighbours == 4)
    {
      std::cerr << "Checking 4 neighbours...\n";
      if(!check_nitrogen_sp3(i, graph.adjacency_list[i], atoms)) { return false; }
    }
    else if(num_neighbours == 3)
    {
      std::cerr << "Checking 3 neighbours...\n";
      if(!check_nitrogen_sp3(i, graph.adjacency_list[i], atoms) &&
         !check_nitrogen_sp2(i, graph.adjacency_list[i], atoms)) { return false; }
    }
    else if(num_neighbours == 2)
    {
      std::cerr << "Checking 2 neighbours...\n";
      if(!check_nitrogen_sp2(i, graph.adjacency_list[i], atoms) &&
         !check_nitrogen_sp(i, graph.adjacency_list[i], atoms)) { return false; }
    }
    else if(num_neighbours == 1) { continue; }
    else { throw std::runtime_error(std::string("STRANGE NITROGEN DETECTED!!  ") + std::to_string(num_neighbours)); }
  }
  return true;
}

bool is_CH3_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph)
{
  if(atoms[target].name != "C") { return false; }
  auto num_H = std::count_if(std::begin(graph.adjacency_list[target]), std::end(graph.adjacency_list[target]),
                             [&atoms](auto&& bond){ return atoms[bond.neighbour].name == "H"; });
  return num_H == 3;
}

bool is_NH2_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph)
{
  if(atoms[target].name != "C") { return false; }
  auto num_H = std::count_if(std::begin(graph.adjacency_list[target]), std::end(graph.adjacency_list[target]),
                             [&atoms](auto&& bond){ return atoms[bond.neighbour].name == "H"; });
  return num_H == 2;
}

bool is_COO_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph)
{
  if(atoms[target].name != "C") { return false; }
  std::size_t single_oxy_neighbours = 0;
  for(auto const& bond : graph.adjacency_list[target])
  {
    // Add if oxygen neighbour has only one neighbour itself
    if(atoms[bond.neighbour].name == "O" and graph.adjacency_list[bond.neighbour].size() == 1) { single_oxy_neighbours += 1; }
  }
  return single_oxy_neighbours == 2;
}

bool is_N_charged_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph)
{
  if(atoms[target].name != "N") { return false; }
  return graph.adjacency_list[target].size() > 3;
}

bool has_metal_bond(molecular_graph const& graph, std::size_t target)
{
  for(auto const& bond : graph.adjacency_list[target])
  {
    if(bond.type == bond_type::METAL) { return true; }
  }
  return false;
}

#ifdef LIBCIF_BUILD_HASH
#include "compute_blake2b_hash.h"
std::string get_label(std::size_t index, std::size_t num_atoms)
{
  auto cell_num = cell_number(index, num_atoms);
  switch (cell_num)
  {
    case 0: return "(0, 0, 0)";

    case 1: return "(-1, -1, -1)"; 
    case 2: return "(-1, -1, 0)"; 
    case 3: return "(-1, -1, 1)"; 
    case 4: return "(-1, 0, -1)"; 
    case 5: return "(-1, 0, 0)"; 
    case 6: return "(-1, 0, 1)"; 
    case 7: return "(-1, 1, -1)"; 
    case 8: return "(-1, 1, 0)"; 
    case 9: return "(-1, 1, 1)"; 

    case 10: return "(0, -1, -1)"; 
    case 11: return "(0, -1, 0)"; 
    case 12: return "(0, -1, 1)"; 
    case 13: return "(0, 0, -1)"; 
    case 14: return "(0, 0, 1)"; 
    case 15: return "(0, 1, -1)"; 
    case 16: return "(0, 1, 0)"; 
    case 17: return "(0, 1, 1)"; 

    case 18: return "(1, -1, -1)"; 
    case 19: return "(1, -1, 0)"; 
    case 20: return "(1, -1, 1)"; 
    case 21: return "(1, 0, -1)"; 
    case 22: return "(1, 0, 0)"; 
    case 23: return "(1, 0, 1)"; 
    case 24: return "(1, 1, -1)"; 
    case 25: return "(1, 1, 0)"; 
    case 26: return "(1, 1, 1)"; 
    default: throw std::runtime_error("nahui");
  }

}

void label_graph(molecular_graph& graph, std::size_t orig_size)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    for(auto& bond : graph.adjacency_list[i])
    {
      bond.label = get_label(bond.neighbour, orig_size);
    }
  }
}

std::string neighbourhood_aggregate(std::size_t target,
                                    molecular_graph const& graph,
                                    std::vector<std::string> const& node_labels,
                                    bool use_edge_attr)
{
  std::vector<std::string> label_list;
  for(auto const& bond : graph.adjacency_list[target])
  {
    std::string prefix = "";
    if(use_edge_attr) { prefix = bond.label; }
    if((!use_edge_attr) || (bond.neighbour>=target))
    {
      label_list.push_back(prefix + node_labels[bond.neighbour]);
    }
  }
  std::sort(std::begin(label_list), std::end(label_list));

  std::string aggregated = node_labels[target];
  for(auto const& label : label_list) { aggregated += label; }
  return aggregated;
}

std::string hash_label(std::string const& str, int digest_size)
{
  std::string output_hex = computeBlake2bHash(str.data(), str.size(), digest_size);
  return output_hex; // return the hex string
}

std::string wl_graph_hash(molecular_graph const& graph,
                          bool use_edge_attr,
                          bool use_node_attr,
                          int num_iteration,
                          std::size_t digest_size)
{
  std::vector<std::string> node_labels(graph.adjacency_list.size());

  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    if(use_node_attr) { node_labels[i] = graph.node_labels[i]; }
    else { node_labels[i] = std::to_string(graph.adjacency_list[i].size()); }
  }

  std::vector<std::pair<std::string, int>> subgraph_hash_counts;
  
  for(int n=0; n<num_iteration; ++n)
  {
    std::vector<std::string> new_labels;
    new_labels.reserve(graph.adjacency_list.size());
    for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
    {
      std::string aggregated = neighbourhood_aggregate(i, graph, node_labels, use_edge_attr);
      new_labels.push_back(hash_label(aggregated, digest_size));
    }

    node_labels = std::move(new_labels);

    std::unordered_map<std::string, int> counts;
    for(auto const& label : node_labels) { counts[label]++; }
    
    std::vector<std::pair<std::string,int>> sorted_counts(counts.begin(), counts.end());
    std::sort(std::begin(sorted_counts), std::end(sorted_counts),
              [](auto&& a, auto&& b){ return a.first < b.first; });

    for(auto const& count : sorted_counts) { subgraph_hash_counts.push_back(count); }
  }

  std::ostringstream oss;
  oss << "(";
  for(std::size_t i=0; i<subgraph_hash_counts.size(); ++i)
  {
    oss << "('" << subgraph_hash_counts[i].first << "', " << subgraph_hash_counts[i].second << ")";
    if(i+1 < subgraph_hash_counts.size()) { oss << ", "; }
  }
  oss << ")";

  std::string hash = hash_label(oss.str(), digest_size);
  return hash;
}
#endif
