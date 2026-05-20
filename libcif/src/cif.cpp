// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "cif.hpp"
#include "linalg.hpp"
#include "utils.hpp"

#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <array>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <iostream>
#include <stdexcept>

double atoms_range(atom const& a, atom const& b)
{
	double r = 0;
	for(int k=0; k<3; ++k)
	{
		r += (a.coord[k]-b.coord[k]) * (a.coord[k]-b.coord[k]);
	}
	return std::sqrt(r);
}

void print_atom(atom const& a, std::ostream* os)
{
  *os << a.name << ' ' << a.coord[0] << ' ' << a.coord[1] << ' ' << a.coord[2] << '\n';
}

void print_xyz(std::vector<atom> const& atoms, std::string const& comment, std::ostream* os)
{
	*os << atoms.size() << '\n';
	*os << comment << '\n';
	for(auto const& a : atoms)
	{
		*os << a.name << ' ' << a.coord[0] << ' ' << a.coord[1] << ' ' << a.coord[2] << '\n';
	}
}

std::pair<std::vector<atom>,std::string> read_xyz(std::istream& is)
{
  std::size_t n;
  is >> n;
  std::string tmp;
  std::string comment_line;
  std::getline(is,tmp);
  std::getline(is, comment_line);
  std::vector<atom> atoms;
  atoms.reserve(n);

  for(std::size_t i=0; i<n; ++i)
  {
    atom a;
    is >> a.name >> a.coord[0] >> a.coord[1] >> a.coord[2];
    atoms.push_back(a);
  }
  return {atoms,comment_line};
}

cif_content::cif_content(std::string const& filename)
{
  auto lines = file_to_lines_strip(filename);
  std::vector<std::string> block_keywords { "data_", "loop_" };
  std::vector<std::size_t> block_starts;
  for(std::size_t i=0; i<lines.size(); ++i)
  {
    for(auto const& block_kw : block_keywords)
    {
      if(lines[i].find(block_kw)!=std::string::npos)
      {
        block_starts.push_back(i);
      }
    }
  }

  for(std::size_t i=0; i<block_starts.size(); ++i)
  {
    cif_block block;
    if(find_kw(lines[block_starts[i]], "data_"))
    {
      block.block_type = cif_block_type::data;
    }
    else if(find_kw(lines[block_starts[i]], "loop_"))
    {
      block.block_type = cif_block_type::loop;
    }

    auto start_block_it = std::next(std::begin(lines), block_starts[i]+1);
    auto end_block_it = std::begin(lines);
    if(i == block_starts.size()-1) { end_block_it = std::end(lines); }
    else
    {
      end_block_it = std::next(std::begin(lines), block_starts[i+1]);
    }
    block.block_content = std::vector<std::string>(start_block_it, end_block_it);

    this->content.push_back(block);
  }

  if(block_starts[0]!=0)
  {
    cif_block implicit_data_block;
    implicit_data_block.block_type = cif_block_type::data;
    implicit_data_block.block_content = std::vector<std::string>(std::begin(lines), std::next(std::begin(lines), block_starts[0]));
    this->content.push_back(implicit_data_block);
  }
}

cif_data::cif_data(std::string const& filename)
{
  auto parsed = cif_content(filename);

  for(auto const& block : parsed.content)
  {
    switch(block.block_type)
    {
      case cif_block_type::data:
        process_data_block(block.block_content);
        break;
      case cif_block_type::loop:
        process_loop_block(block.block_content);
        break;
    }
  }

  if(this->num_atoms == 0) { throw std::runtime_error("ERROR! No atoms found!"); }

  for(auto& atom : fract_atoms)
  {
    while(atom.coord[0] < 0) { atom.coord[0] += 1; }
    while(atom.coord[1] < 0) { atom.coord[1] += 1; }
    while(atom.coord[2] < 0) { atom.coord[2] += 1; }
    while(atom.coord[0] >= 1) { atom.coord[0] -= 1; }
    while(atom.coord[1] >= 1) { atom.coord[1] -= 1; }
    while(atom.coord[2] >= 1) { atom.coord[2] -= 1; }
  }
}

void cif_data::process_data_block(std::vector<std::string> const& data)
{
  std::vector<std::string> const data_keywords
  {
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
    "_symmetry_space_group_name_H-M",
    "_chemical_formula_sum"
  };

  std::unordered_map<std::string, std::string> data_map {};
  for(auto const& line : data)
  {
    std::string key{};
    std::string value{};
    std::istringstream iss(line);
    iss >> key;
    if(next_non_blank_in_stream(iss) == '\'') { iss >> std::quoted(value, '\''); }
    else if(next_non_blank_in_stream(iss) == '\"') { iss >> std::quoted(value, '\"'); }
    else { iss >> value; }
    data_map[key] = value;
  }

  // Note: only checks within a single block — required keywords spread across
  // multiple blocks will not be merged.
  for(std::size_t i=0; i<6; ++i)
  {
    if(data_map.count(data_keywords[i]) == 0) { return; }
  }

  try
  {
    double a = std::stod(data_map["_cell_length_a"]);
    double b = std::stod(data_map["_cell_length_b"]);
    double c = std::stod(data_map["_cell_length_c"]);
    double alpha = std::stod(data_map["_cell_angle_alpha"]);
    double beta = std::stod(data_map["_cell_angle_beta"]);
    double gamma = std::stod(data_map["_cell_angle_gamma"]);
    this->cell = make_unit_cell(a, b, c, alpha, beta, gamma);
  }
  catch (std::invalid_argument const&)
  {
    std::cerr << "ERROR: Can't read data block: wrong number format!" << std::endl;
    std::cout << "BAD" << std::endl;
    std::exit(1);
  }
  if(data_map.contains("_chemical_formula_sum")) { this->formula = data_map["_chemical_formula_sum"]; }
}

void cif_data::process_loop_block(std::vector<std::string> const& data)
{
  std::vector<std::string> const atom_keywords
  {
    "_atom_site_type_symbol",
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
    "_atom_site_Cartn_x",
    "_atom_site_Cartn_y",
    "_atom_site_Cartn_z",
    "_atom_site_charge"
  };


  auto end_kws_it = std::find_if(std::begin(data), std::end(data),
                                 [](auto const& l){ return first_non_blank(l) != '_'; });
  std::size_t num_kws = std::distance(std::begin(data), end_kws_it);

  // Check if this loop block contains other data block
  auto end_table_it = std::find_if(end_kws_it, std::end(data),
                                 [](auto const& l){ return first_non_blank(l) == '_'; });
  if(end_table_it != std::end(data))
  {
    process_data_block(std::vector<std::string>(end_table_it, std::end(data)));
    process_loop_block(std::vector<std::string>(std::begin(data), end_table_it));
    return;
  }



  std::vector<std::string> kws;
  for(std::size_t i=0; i<num_kws; ++i)
  {
    std::istringstream iss(data[i]);
    std::string tmp;
    iss >> tmp;
    kws.push_back(tmp);
  }

  std::unordered_map<std::string, std::vector<std::string>> loop_map;

  for(std::size_t i=num_kws; i<data.size(); ++i)
  {
    std::vector<std::string> tmp_line;
    std::string tmp;
    std::istringstream iss(data[i]);
    for(std::size_t j=0; j<num_kws; ++j)
    {
      iss >> tmp;
      tmp_line.push_back(tmp);
    }

    for(std::size_t j=0; j<num_kws; ++j)
    {
      loop_map[kws[j]].push_back(tmp_line[j]);
    }
  }

  if(loop_map.count(atom_keywords[0]) == 0) { return; }
  this->num_atoms = loop_map[atom_keywords[0]].size();

  if(loop_map.count(atom_keywords[1]))
  {
    this->fract_atoms.resize(num_atoms);
    for(std::size_t i=0; i<num_atoms; ++i)
    {
      this->fract_atoms[i].coord[0] = std::stod(loop_map[atom_keywords[1]][i]);
      this->fract_atoms[i].coord[1] = std::stod(loop_map[atom_keywords[2]][i]);
      this->fract_atoms[i].coord[2] = std::stod(loop_map[atom_keywords[3]][i]);
    }
  }
  
  if(loop_map.count(atom_keywords[4]))
  {
    this->abs_atoms.resize(num_atoms);
    for(std::size_t i=0; i<num_atoms; ++i)
    {
      this->abs_atoms[i].coord[0] = std::stod(loop_map[atom_keywords[4]][i]);
      this->abs_atoms[i].coord[1] = std::stod(loop_map[atom_keywords[5]][i]);
      this->abs_atoms[i].coord[2] = std::stod(loop_map[atom_keywords[6]][i]);
    }
  }

  if(this->fract_atoms.size() != 0)
  {
    for(std::size_t i=0; i<num_atoms; ++i) { this->fract_atoms[i].name = loop_map[atom_keywords[0]][i]; }
  }
  if(this->abs_atoms.size() != 0)
  {
    for(std::size_t i=0; i<num_atoms; ++i) { this->abs_atoms[i].name = loop_map[atom_keywords[0]][i]; }
  }

  if(this->fract_atoms.size() != 0)
  {
    if(loop_map.count(atom_keywords[7]))
    {
      for(std::size_t i=0; i<num_atoms; ++i) { this-> fract_atoms[i].q = std::stod(loop_map[atom_keywords[7]][i]); }
    }
  }
  if(this->abs_atoms.size() != 0)
  {
    if(loop_map.count(atom_keywords[7]))
    {
      for(std::size_t i=0; i<num_atoms; ++i) { this-> abs_atoms[i].q = std::stod(loop_map[atom_keywords[7]][i]); }
    }
  }


}

std::vector<atom> const& cif_data::get_any_atoms_ref() const
{
  if(fract_atoms.size() != 0) { return std::cref(fract_atoms); }
  else if(abs_atoms.size() != 0) { return std::cref(abs_atoms); }
  else { throw std::runtime_error("No atoms to return!!!!"); }
}

void cif_data::set_unit_cell()
{
  this->cell = unit_cell();
}

atom fract_to_abs(atom const& fract_atom, unit_cell const& cell)
{
  matrix33 m;
  m(0,0) = cell.axes.a[0];
  m(1,0) = cell.axes.a[1];
  m(2,0) = cell.axes.a[2];
  
  m(0,1) = cell.axes.b[0];
  m(1,1) = cell.axes.b[1];
  m(2,1) = cell.axes.b[2];
  
  m(0,2) = cell.axes.c[0];
  m(1,2) = cell.axes.c[1];
  m(2,2) = cell.axes.c[2];

  atom abs_atom;
  abs_atom.name = fract_atom.name;
  abs_atom.coord = multiply(m, fract_atom.coord);

  return abs_atom;
}

atom abs_to_fract(atom const& abs_atom, unit_cell const& cell)
{
  matrix33 m;
  m(0,0) = cell.axes.a[0];
  m(1,0) = cell.axes.a[1];
  m(2,0) = cell.axes.a[2];
  
  m(0,1) = cell.axes.b[0];
  m(1,1) = cell.axes.b[1];
  m(2,1) = cell.axes.b[2];
  
  m(0,2) = cell.axes.c[0];
  m(1,2) = cell.axes.c[1];
  m(2,2) = cell.axes.c[2];

  auto m_inv = inverse(m);

  atom fract_atom;
  fract_atom.name = abs_atom.name;
  fract_atom.coord = multiply(m_inv, abs_atom.coord);

  return fract_atom;
}

void cif_data::calc_abs_coords()
{
  // Don't calculate if already exists
  if(this->abs_atoms.size() != 0) { return; }

  this->abs_atoms = this->fract_atoms;

  for(std::size_t i=0; i<fract_atoms.size(); ++i)
  {
    abs_atoms[i].coord[0] = fract_atoms[i].coord[0]*this->cell.axes.a[0] +
                            fract_atoms[i].coord[1]*this->cell.axes.b[0] +
                            fract_atoms[i].coord[2]*this->cell.axes.c[0];
    abs_atoms[i].coord[1] = fract_atoms[i].coord[0]*this->cell.axes.a[1] +
                            fract_atoms[i].coord[1]*this->cell.axes.b[1] +
                            fract_atoms[i].coord[2]*this->cell.axes.c[1];
    abs_atoms[i].coord[2] = fract_atoms[i].coord[0]*this->cell.axes.a[2] +
                            fract_atoms[i].coord[1]*this->cell.axes.b[2] +
                            fract_atoms[i].coord[2]*this->cell.axes.c[2];
  }
}

void cif_data::calc_fract_coords()
{
  if(this->fract_atoms.size() != 0) { return; }

  this->fract_atoms = this->abs_atoms;

  for(std::size_t i=0; i<abs_atoms.size(); ++i)
  {
    this->fract_atoms[i] = abs_to_fract(abs_atoms[i], this->cell);
  }
}


void cif_data::calc_supercell(int x_min, int x_max,
                              int y_min, int y_max,
                              int z_min, int z_max,
                              bool calc_fract)
{
  // Don't calculate if already exists
  if(this->supercell.size() != 0) { return; }
	this->calc_abs_coords();

	std::size_t const n_atoms = this->abs_atoms.size();
  this->supercell = this->abs_atoms;
  int const n_cells = (x_max-x_min+1) * (y_max-y_min+1) * (z_max-z_min+1);
	this->supercell.reserve(n_atoms*n_cells);

	for(int i=x_min; i<=x_max; ++i)
	{
		for(int j=y_min; j<=y_max; ++j)
		{
			for(int k=z_min; k<=z_max; ++k)
			{
				if(i==0 and j==0 and k==0) { continue; }

				std::vector<atom> tmp(this->abs_atoms);
				for(std::size_t n=0; n<n_atoms; ++n)
				{
					tmp[n].coord[0] += i * this->cell.axes.a[0] +
                             j * this->cell.axes.b[0] +
                             k * this->cell.axes.c[0];

					tmp[n].coord[1] += i * this->cell.axes.a[1] +
                             j * this->cell.axes.b[1] +
                             k * this->cell.axes.c[1];

					tmp[n].coord[2] += i * this->cell.axes.a[2] +
                             j * this->cell.axes.b[2] +
                             k * this->cell.axes.c[2];
				}
        append_vector(tmp, this->supercell);
			}
		}
	}
  if(!calc_fract) { return; }
  this->fract_supercell = this->fract_atoms;
  this->fract_supercell.reserve(n_atoms*n_cells);
  for(int i=x_min; i<=x_max; ++i)
  {
		for(int j=y_min; j<=y_max; ++j)
		{
			for(int k=z_min; k<=z_max; ++k)
			{
				if(i==0 and j==0 and k==0) { continue; }

        std::vector<atom> tmp(this->fract_atoms);
        for(std::size_t n=0; n<n_atoms; ++n)
        {
          tmp[n].coord[0] += i;
          tmp[n].coord[1] += j;
          tmp[n].coord[2] += k;
        }
        append_vector(tmp, this->fract_supercell);
      }
    }
  }
}

void cif_data::dump_primcell(std::ostream* os) const
{
  *os << this->cell.a << ' ' << this->cell.b << ' ' << this->cell.c << '\n';
  *os << this->cell.alpha << ' ' << this->cell.beta << ' ' << this->cell.gamma << '\n';
  print_xyz(this->abs_atoms, "comment", os);
}

void chromosome::print_mutations(std::ostream* os) const
{
  for(auto const& g : this->genome)
  {
    *os << g.mutation_number << ' ';
  }
  *os << '\n';
}

void chromosome::print_full(std::ostream* os) const
{
  print_mutations(os);

  for(auto const& gene : genome)
  {
    *os << gene.mutation_number << ' ' << gene.sites.size() << '\n';
    for(auto const& sites : gene.sites)
    {
      for(auto const& idx : sites)
      {
        *os << idx << ' ';
      }
      *os << '\n';
    }
  }
}

bool is_inside_cell(atom const& abs_atom, unit_cell const& cell)
{
  auto fract_atom = abs_to_fract(abs_atom, cell);

  bool a_inside = ((fract_atom.coord[0] >=0) and (fract_atom.coord[0] <=1));
  bool b_inside = ((fract_atom.coord[1] >=0) and (fract_atom.coord[1] <=1));
  bool c_inside = ((fract_atom.coord[2] >=0) and (fract_atom.coord[2] <=1));

  return (a_inside and b_inside and c_inside);
}

void print_graph(molecular_graph const& graph, std::ostream* os)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    *os << i + 1 << '\t';
    for(auto const& bond : graph.adjacency_list[i])
    {
        *os << bond.neighbour + 1 << ' ';
    }
    *os << '\n';
  }
}

void print_graph_v2(molecular_graph const& graph, std::ostream* os)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    for(auto const& bond : graph.adjacency_list[i])
    {
      *os << i << ' ' << bond.neighbour << ' ' << bond.label << std::endl;
    }
  }
}

void print_graph_sane(molecular_graph const& graph, std::ostream* os)
{
  for(std::size_t i=0; i<graph.adjacency_list.size(); ++i)
  {
    for(auto const& bond : graph.adjacency_list[i])
    {
      *os << i << ' ' << bond.neighbour << ' ' << bond.r << ' '
        << bond.label << ' ' << static_cast<int>(bond.type) << '\n';
    }
  }
}

void print_cif(cif_data const& cif, std::ostream* out)
{
  *out << "data_\n";
  *out << "_cell_length_a " << cif.cell.a << '\n';
  *out << "_cell_length_b " << cif.cell.b << '\n';
  *out << "_cell_length_c " << cif.cell.c << '\n';
  *out << "_cell_angle_alpha " << cif.cell.alpha << '\n';
  *out << "_cell_angle_beta " << cif.cell.beta << '\n';
  *out << "_cell_angle_gamma " << cif.cell.gamma << '\n';
  *out << "loop_\n";
  if(cif.fract_atoms.size() != 0)
  {
    *out << "_atom_site_type_symbol\n";
    *out << "_atom_site_fract_x\n";
    *out << "_atom_site_fract_y\n";
    *out << "_atom_site_fract_z\n";
    for(auto const& atom : cif.fract_atoms)
    {
      *out << atom.name << ' ' << atom.coord[0] << ' ' <<
        atom.coord[1] << ' ' << atom.coord[2] << '\n';
    }
  }
  else if(cif.abs_atoms.size() != 0)
  {
    *out << "_atom_site_type_symbol\n";
    *out << "_atom_site_Cartn_x\n";
    *out << "_atom_site_Cartn_y\n";
    *out << "_atom_site_Cartn_z\n";
    for(auto const& atom : cif.abs_atoms)
    {
      *out << atom.name << ' ' << atom.coord[0] << ' ' <<
        atom.coord[1] << ' ' << atom.coord[2] << '\n';
    }
  }
  else { throw std::runtime_error("ERROR! No atoms to print!"); }
}

void reflect_if_needed(atom& abs_atom, unit_cell const& cell)
{
  auto fract_atom = abs_to_fract(abs_atom, cell);

  if(fract_atom.coord[0] > 1.0) { fract_atom.coord[0] -= 1; }
  if(fract_atom.coord[0] < 0.0) { fract_atom.coord[0] += 1; }

  if(fract_atom.coord[1] > 1.0) { fract_atom.coord[1] -= 1; }
  if(fract_atom.coord[1] < 0.0) { fract_atom.coord[1] += 1; }

  if(fract_atom.coord[2] > 1.0) { fract_atom.coord[2] -= 1; }
  if(fract_atom.coord[2] < 0.0) { fract_atom.coord[2] += 1; }

  abs_atom = fract_to_abs(fract_atom, cell);
}

std::array<std::array<double,3>,3> make_axes_matrix(crystal_axes const& axes)
{
  std::array<std::array<double,3>,3> matrix;
  matrix[0] = axes.a;
  matrix[1] = axes.b;
  matrix[2] = axes.c;
  return matrix;
}
