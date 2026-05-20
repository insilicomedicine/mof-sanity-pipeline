// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"

#include <string>
#include <unordered_set>
#include <utility>
#include <unordered_map>

struct atom_with_radius : public atom
{
	double vdw_radius;

  atom_with_radius(atom a, double r) : atom(std::move(a)), vdw_radius(r) {};
};

struct coordination_spheres
{
	std::string name;
	std::size_t n02;
	std::size_t n04;
	std::size_t n06;
	std::size_t n08;
	std::size_t n10;
};

// Operations on cifs
std::string calc_cif_formula(cif_data const& cif);



bool cif_is_good(cif_data& cif,
                 double r_min = 0.1, double r_max = 5.0,
                 double tolerance = 1.0,
                 std::string const& strategy = "cordero");

bool cif_is_good_v2(cif_data& cif,
                 double r_min = 0.1, double r_max = 5.0,
                 double tolerance = 1.0,
                 std::string const& strategy = "cordero");

double cif_quality(cif_data& cif,
                   double r_min = 0.1, double r_max = 5.0,
                   double tolerance = 1.0,
                   std::string const& strategy = "cordero",
                   double penalty = 0.5);

std::string calc_cif_hash(cif_data& cif,
                          bool use_edge_attr, bool use_node_attr,
                          double tolerance = 1.0,
                          std::string const& strategy = "cordero");

double calc_PSA_ratio(cif_data const& cif, double thr, std::string const& method);

std::vector<linker> separate_linkers(cif_data& cif,
                                     double r_min = 0.1, double r_max = 5.0,
                                     double tolerance = 1.0,
                                     std::string const& strategy = "cordero");

std::vector<linker> separate_linkers_v2(cif_data& cif, std::size_t min_size,
                                        double r_min = 0.1, double r_max = 5.0,
                                        double tolerance = 1.0,
                                        std::string const& strategy = "cordero");

std::unordered_set<std::string> separate_linkers_smiles(cif_data& cif,
                                                        double r_min = 0.1, double r_max = 5.0,
                                                        double tolerance = 1.0,
                                                        std::string const& strategy = "cordero");
std::unordered_set<std::string> separate_linkers_smiles_v2(cif_data& cif, std::size_t min_size,
                                                        double r_min = 0.1, double r_max = 5.0,
                                                        double tolerance = 1.0,
                                                        std::string const& strategy = "cordero");

void decompose_mof(cif_data& cif,
                   double tolerance = 1.3,
                   std::string const& strategy = "cordero",
                   bool use_linear_graph = false);

std::vector<coordination_spheres> calc_metals_coord_num(cif_data const& cif);

double calc_volume(cif_data const& cif);
double calc_mass(cif_data const& cif);
double calc_density(cif_data const& cif);
double calc_total_charge(cif_data const& cif);

std::vector<atom> indices_to_atoms(std::vector<std::size_t> const& indices, std::vector<atom> const& orig_atoms);

bool all_atoms_unsaturated(std::vector<std::size_t> const& targets,
                           std::vector<atom> const& atoms, molecular_graph const& graph);

chromosome make_chromosome(cif_data& cif, double r_min, double r_max, double tolerance, std::string const& strategy);

void print_mutant(cif_data const& cif, std::ostream* os = &std::cout);
cif_data read_mutant(std::istream& is);

void mutate(cif_data& cif, std::vector<int> const& new_mutations);

std::tuple<std::vector<std::string>,int> calc_metal_charges(cif_data& cif,
                                                            double r_min = 0.1, double r_max = 5.0,
                                                            double tolerance = 1.0,
                                                            std::string const& strategy = "cordero");

std::unordered_map<std::string, std::vector<linker>> group_by_smiles(std::vector<linker> const& linkers);

void cell_prim_spg(cif_data& cif, double tolerance = 0.25);

void remove_all_outside(std::vector<linker>& structs, std::size_t orig_size);
