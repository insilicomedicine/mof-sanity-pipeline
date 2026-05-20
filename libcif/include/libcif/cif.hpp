// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "linalg.hpp"
#include <array>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

struct atom
{
  vec3d coord;
  double q = 0; // charge
  std::string name;
};

double atoms_range(atom const& a, atom const& b);

void print_atom(atom const& a, std::ostream* os = &std::cout);
void print_xyz(std::vector<atom> const& atoms, std::string const& comment = "Comment line", std::ostream* os = &std::cout);
std::pair<std::vector<atom>, std::string> read_xyz(std::istream& is);

enum class cif_block_type { data, loop };

struct cif_block
{
  cif_block_type block_type;
  std::vector<std::string> block_content;
};

struct cif_content
{
  std::vector<cif_block> content;

  cif_content(std::string const& filename);
};

enum class bond_type { SINGLE=1, DOUBLE=2, TRIPLE=3, METAL=-1 };

struct bond
{
  std::size_t neighbour;
  double r;
  bond_type type;
  std::string label = "";
};

struct molecular_graph
{
  std::vector<std::vector<bond>> adjacency_list;
  std::vector<std::string> node_labels;
  std::string hash;
};

struct graph_params
{
  double r_min = 0.5;
  double r_max = 5.0;
  double tolerance = 1.2;
};

struct gene
{
  int mutation_number;
  std::vector<std::vector<std::size_t>> sites;
};

struct chromosome
{
  std::vector<gene> genome;

  void print_mutations(std::ostream* os = &std::cout) const;
  void print_full(std::ostream* os = &std::cout) const;
};

struct crystal_axes
{
  std::array<double,3> a;
  std::array<double,3> b;
  std::array<double,3> c;
};

std::array<std::array<double,3>,3> make_axes_matrix(crystal_axes const& axes);

struct unit_cell
{
  double a;
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;

  crystal_axes axes;
};

constexpr unit_cell make_unit_cell(double a, double b, double c, double alpha, double beta, double gamma)
{
	double const alpha_rad = alpha * M_PI / 180;
	double const beta_rad = beta * M_PI / 180;
	double const gamma_rad = gamma * M_PI / 180;

  crystal_axes axes;

	axes.a[0] = a;
	axes.a[1] = 0;
	axes.a[2] = 0;
	axes.b[0] = b * cos(gamma_rad);
	axes.b[1] = b * sin(gamma_rad);
	axes.b[2] = 0;
	axes.c[0] = c * cos(beta_rad);
	axes.c[1] = c * ((cos(alpha_rad) - cos(gamma_rad)*cos(beta_rad)) / sin(gamma_rad));
	axes.c[2] = sqrt(c*c - axes.c[0]*axes.c[0] - axes.c[1]*axes.c[1]);

  return {a,b,c,alpha,beta,gamma,axes};
}

inline unit_cell make_unit_cell(crystal_axes const& axes)
{
  auto a_vec = vec3d(axes.a);
  auto b_vec = vec3d(axes.b);
  auto c_vec = vec3d(axes.c);

  double a = norm(a_vec);
  double b = norm(b_vec);
  double c = norm(c_vec);

  double ab = scalar_product(a_vec, b_vec);
  double ac = scalar_product(a_vec, c_vec);
  double bc = scalar_product(b_vec, c_vec);

  double alpha_rad = std::acos(bc / (b*c));
  double beta_rad = std::acos(ac/(a*c));
  double gamma_rad = std::acos(ab/(a*b));

  double const alpha = alpha_rad * 180.0 / M_PI;
  double const beta = beta_rad * 180.0 / M_PI;
  double const gamma = gamma_rad * 180.0 / M_PI;

  return {a, b, c, alpha, beta, gamma, axes};
}

struct cif_data
{
  unit_cell cell;

  std::size_t num_atoms;
	std::string formula; // Brute formula

	std::vector<atom> fract_atoms;
  std::vector<atom> abs_atoms;
  std::vector<atom> supercell;
  std::vector<atom> fract_supercell;

  molecular_graph graph;

  chromosome cmsm;

  cif_data() = default;
	cif_data(std::string const& filename);

  void set_unit_cell();

  void calc_abs_coords();
  void calc_fract_coords();
  void calc_supercell(int x_min, int x_max,
                      int y_min, int y_max,
                      int z_min, int z_max,
                      bool calc_fract = true);

  std::vector<atom> const& get_any_atoms_ref() const;

  void dump_primcell(std::ostream* os = &std::cout) const;

private:
	void process_data_block(std::vector<std::string> const& data);
	void process_loop_block(std::vector<std::string> const& data);
};


struct linker
{
  std::vector<std::size_t> indices = {};
  std::vector<atom> atoms = {};
  std::string smiles = {};
  molecular_graph graph = {};
  std::vector<std::size_t> canon_order = {};
  std::vector<std::size_t> orig_indices = {};
  std::string comment = {};
  bool infinite = false;
};


bool is_inside_cell(atom const& abs_atom, unit_cell const& cell);

std::array<double,3> determine_cell_translation(std::size_t n, std::size_t n_atoms);

void print_cif(cif_data const& cif, std::ostream* out = &std::cout);
void print_graph(molecular_graph const& graph, std::ostream* os = &std::cout);
void print_graph_v2(molecular_graph const& graph, std::ostream* os);
void print_graph_sane(molecular_graph const& graph, std::ostream* os);

void reflect_if_needed(atom& abs_atom, unit_cell const& cell);
