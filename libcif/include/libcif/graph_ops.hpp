// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"
#include <unordered_set>


std::vector<bond> find_bonds(std::size_t target, std::vector<atom> const& atoms,
                                  double tolerance, std::string const& strategy);

molecular_graph calc_connect_graph(std::vector<atom> const& atoms, std::size_t graph_size,
                                   double tolerance, std::string const& strategy);

molecular_graph calc_connect_graph_v2(std::vector<atom> const& atoms, std::size_t graph_size,
                                      double tolerance, std::string const& strategy);

molecular_graph calc_connect_graph_v3(std::vector<atom> const& atoms, std::vector<atom> const& fract_atoms,
                                      std::size_t graph_size, unit_cell const& cell,
                                      double tolerance, std::string const& strategy);

bool check_overlap(molecular_graph const& graph, std::vector<atom> const& atoms, double r_min = 0.5);

void remove_too_long_bonds(molecular_graph& graph, double r_max = 5.0);

bool calc_graph(cif_data& cif,
                double r_min, double r_max,
                double tolerance,
                std::string const& strategy);

bool verify_trim_graph(molecular_graph& graph, std::vector<atom> const& atoms, double r_min = 0.5, double r_max = 5.0);

bool all_atoms_unsaturated(std::vector<std::size_t> const& targets,
                           std::vector<atom> const& atoms, molecular_graph const& graph);

std::vector<std::vector<std::size_t>> find_all_n_member_rings(molecular_graph const& graph, std::size_t n);

void make_linker_dfs_v2(std::size_t target,
                        molecular_graph const& graph,
                        std::vector<atom> const& atoms,
                        std::vector<std::size_t>& linker,
                        std::unordered_set<std::size_t>& visited);

void make_node_dfs_v2(std::size_t target,
                      molecular_graph const& graph,
                      std::vector<atom> const& atoms,
                      std::vector<std::size_t>& node,
                      std::unordered_set<std::size_t>& visited,
                      std::unordered_set<std::size_t> const& linkers_indices);

void make_disconnected_dfs_v2(std::size_t target,
                      molecular_graph const& graph,
                      std::vector<atom> const& atoms,
                      std::vector<std::size_t>& node,
                      std::unordered_set<std::size_t>& visited,
                      std::unordered_set<std::size_t> const& linkers_indices,
                      std::unordered_set<std::size_t> const& nodes_indices);

std::vector<std::size_t> find_overlap_atoms(molecular_graph const& graph, double r_min = 0.5);

void find_connected_dfs(std::size_t target, std::unordered_set<std::size_t>& visited, std::vector<std::size_t>& subgraph, molecular_graph const& graph);

molecular_graph transform_to_original_indices(molecular_graph graph, std::size_t num_atoms);

bool check_disconnected(cif_data const& cif, std::size_t max_err_size = 3);

std::vector<std::vector<std::size_t>> find_subgraphs(molecular_graph const& graph);

bool check_carbons(molecular_graph const& graph, std::vector<atom> const& atoms);
bool check_nitrogens(molecular_graph const& graph, std::vector<atom> const& atoms);

void label_graph(molecular_graph& graph, std::size_t orig_size);

bool has_metal_bond(molecular_graph const& graph, std::size_t target);

bool is_CH3_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph);
bool is_NH2_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph);
bool is_COO_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph);
bool is_N_charged_group(std::size_t target, std::vector<atom> const& atoms, molecular_graph const& graph);

std::string wl_graph_hash(molecular_graph const& graph,
                          bool use_edge_attr = false,
                          bool use_node_attr = false,
                          int num_iteration = 6,
                          std::size_t digest_size = 16);
