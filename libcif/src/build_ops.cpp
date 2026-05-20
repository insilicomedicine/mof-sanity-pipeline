// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "build_ops.hpp"
#include "building_blocks.hpp"
#include "linalg.hpp"
#include "cif.hpp"
#include "utils.hpp"
#include <sstream>
#include <stdexcept>
#include <unordered_set>

[[nodiscard]] building_blocks::edge_type read_xyz_to_edge(std::string const& filename)
{
  std::fstream inp(filename);
  if(inp.fail()) { throw std::runtime_error("Can't open edge file!"); }
  auto [atoms,comment] = read_xyz(inp);

  std::istringstream iss(comment);
  std::array<double,3> normal;
  std::vector<int> connection_points;

  iss >> normal[0] >> normal[1] >> normal[2];

  int tmp;
  while(iss>>tmp)
  {
    connection_points.push_back(tmp);
  }

  return {.atoms=atoms, .connection_points=connection_points, .normal=normal};
}

[[nodiscard]] building_blocks::node_type read_xyz_to_node(std::string const& filename)
{
  std::fstream inp(filename);
  if(inp.fail()) { throw std::runtime_error("Can't open edge file!"); }
  auto [atoms,comment] = read_xyz(inp);

  std::vector<atom> ghosts;
  auto new_end = std::find_if(std::begin(atoms), std::end(atoms), [](atom const& a){ return a.name == "X"; });
  std::copy(new_end, std::end(atoms), std::back_inserter(ghosts));
  atoms.erase(new_end, std::end(atoms));

  return {.atoms=atoms, .ghosts=ghosts};
}

[[nodiscard]] std::vector<building_blocks::node_type> generate_nodes(crystal_axes const& axes, building_blocks::topology const& topology)
{
  std::vector<building_blocks::node_type> nodes;
  for(auto const& node : topology.nodes)
  {
    auto blueprint = node.get_blueprint();
    rotate_around_origin(blueprint.atoms, node.orientation);
    rotate_around_origin(blueprint.ghosts, node.orientation);

    auto new_position = node.frac_coord[0]*vec3d{axes.a} +
                                  node.frac_coord[1]*vec3d{axes.b} +
                                  node.frac_coord[2]*vec3d{axes.c};

    auto translation = new_position - vec3d{blueprint.center.coord};

    translate(blueprint.atoms, translation);
    translate(blueprint.ghosts, translation);

    nodes.push_back(blueprint);
  }

  return nodes;
}

void move_node(building_blocks::node_type& node, building_blocks::edge_type const& edge,
               std::size_t ghost, std::size_t link)
{
  vec3d initial_point = node.ghosts[ghost].coord;
  vec3d final_point = edge.atoms[link].coord;
  auto translation = final_point - initial_point;
  translate(node.atoms, translation);
  translate(node.ghosts, translation);
}

void move_edge(building_blocks::edge_type& edge, building_blocks::node_type const& node,
               std::size_t link, std::size_t ghost)
{
  vec3d initial_point = edge.atoms[link].coord;
  vec3d final_point = node.ghosts[ghost].coord;
  auto translation = final_point - initial_point;
  translate(edge.atoms, translation);
}

[[nodiscard]] cif_data build_hum(building_blocks::topology const& topology)
{
  cif_data cif;
  cif.cell = make_unit_cell(100, 100, 100, 90, 90, 90);

  std::vector<building_blocks::node_type> real_nodes = generate_nodes(cif.cell.axes, topology);

  std::vector<building_blocks::edge_type> real_edges;
  
  for(auto const& edge : topology.edges)
  {
    auto blueprint = edge.get_blueprint();
    if(edge.connections.size() != 2) { throw std::runtime_error("Non-linear linkers not implemented!!"); }
    vec3d orientation_to = vec3d(real_nodes[edge.connections[1].node_number].ghosts[edge.connections[1].ghost_number].coord) -
                           vec3d(real_nodes[edge.connections[0].node_number].ghosts[edge.connections[0].ghost_number].coord);
    vec3d orientation_from = vec3d(blueprint.atoms[blueprint.connection_points[1]].coord) -
                             vec3d(blueprint.atoms[blueprint.connection_points[0]].coord);
    vec3d normal_to = vec3d(edge.normal_to);
    vec3d normal_from = vec3d(blueprint.normal);

    auto rotation_angles = calc_euler_angles(orientation_from, normal_from,
                                             orientation_to, normal_to);
    rotate_around_origin(blueprint.atoms, rotation_angles);

    vec3d translation = vec3d(real_nodes[edge.connections[0].node_number].ghosts[edge.connections[0].ghost_number].coord) -
                        vec3d(blueprint.atoms[blueprint.connection_points[0]].coord);
    translate(blueprint.atoms, translation);
    real_edges.push_back(blueprint);
  }

  for(std::size_t i=1; i<real_nodes.size(); ++i)
  {
    vec3d initial_point = vec3d(real_nodes[i].ghosts[topology.edges[i-1].connections[1].ghost_number].coord);
    vec3d final_point = vec3d(real_edges[i-1].atoms[real_edges[i-1].connection_points[1]].coord);
    vec3d translation = final_point - initial_point;
    translate(real_nodes[i].atoms, translation);
  }

  for(auto const& node : real_nodes)
  {
    append_vector(node.atoms, cif.abs_atoms);
  }
  for(auto const& edge : real_edges)
  {
    append_vector(edge.atoms, cif.abs_atoms);
  }

  double new_a = atoms_range(real_nodes[0].atoms[0], real_nodes[1].atoms[0]);
  double new_b = atoms_range(real_nodes[0].atoms[0], real_nodes[3].atoms[0]);
  double new_c = atoms_range(real_nodes[0].atoms[0], real_nodes[5].atoms[0]);

  cif.cell.a = new_a;
  cif.cell.b = new_b;
  cif.cell.c = new_c;

  return cif;
}

[[nodiscard]] cif_data build_999(building_blocks::topology const& topology, double const c)
{
  cif_data cif;
  cif.cell = make_unit_cell(100, 100, c, 90, 90, 60);

  std::vector<building_blocks::node_type> real_nodes = generate_nodes(cif.cell.axes, topology);
  std::vector<building_blocks::edge_type> real_edges;

  for(auto const& edge : topology.edges)
  {
    auto blueprint = edge.get_blueprint();
    if(edge.connections.size() != 2) { throw std::runtime_error("Non-linear linkers not implemented!!"); }
    vec3d orientation_to = vec3d(real_nodes[edge.connections[1].node_number].ghosts[edge.connections[1].ghost_number].coord) -
                           vec3d(real_nodes[edge.connections[0].node_number].ghosts[edge.connections[0].ghost_number].coord);
    vec3d orientation_from = vec3d(blueprint.atoms[blueprint.connection_points[1]].coord) -
                             vec3d(blueprint.atoms[blueprint.connection_points[0]].coord);
    vec3d normal_to = vec3d(edge.normal_to);
    vec3d normal_from = vec3d(blueprint.normal);

    auto rotation_angles = calc_euler_angles(orientation_from, normal_from,
                                             orientation_to, normal_to);
    rotate_around_origin(blueprint.atoms, rotation_angles);
    real_edges.push_back(blueprint);
  }

  std::unordered_set<std::size_t> visited;
  visited.insert(0);

  for(std::size_t i=0; i<topology.edges.size(); ++i)
  {
    std::size_t node_to_move = topology.edges[i].connections[0].node_number;
    std::size_t target_ghost = topology.edges[i].connections[0].ghost_number;
    if(!visited.contains(node_to_move))
    {
      move_node(real_nodes[node_to_move], real_edges[i], target_ghost, real_edges[i].connection_points[0]);
      visited.insert(node_to_move);
    }

    move_edge(real_edges[i], real_nodes[node_to_move], real_edges[i].connection_points[0], target_ghost);
    
    node_to_move = topology.edges[i].connections[1].node_number;
    target_ghost = topology.edges[i].connections[1].ghost_number;
    if(!visited.contains(node_to_move))
    {
      move_node(real_nodes[node_to_move], real_edges[i], target_ghost, real_edges[i].connection_points[1]);
      visited.insert(node_to_move);
    }
  }

  double new_a = atoms_range(real_nodes[0].atoms[0], real_nodes[4].atoms[0]);
  double new_b = atoms_range(real_nodes[0].atoms[0], real_nodes[5].atoms[0]);

  cif.cell.a = new_a;
  cif.cell.b = new_b;

  append_vector(real_nodes[0].atoms, cif.abs_atoms);
  append_vector(real_nodes[1].atoms, cif.abs_atoms);

  append_vector(real_edges[0].atoms, cif.abs_atoms);
  append_vector(real_edges[3].atoms, cif.abs_atoms);
  append_vector(real_edges[4].atoms, cif.abs_atoms);

  return cif;
}
