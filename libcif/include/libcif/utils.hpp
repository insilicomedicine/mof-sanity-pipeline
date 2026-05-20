// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <vector>

constexpr double const amu_to_g = 1.6605391e-24;
constexpr double const angstrom3_to_cm3 = 1e-24;

inline char first_non_blank(std::string const& str)
{
	return *std::find_if(std::begin(str), std::end(str), [](char c){ return !std::isspace(c); });
}

inline bool find_kw(std::string const& line, std::string const& kw)
{
  if(line.size()==0) { return false; }
  if(line.find(kw)!=std::string::npos) { return true; }
  return false;
};

inline char next_non_blank_in_stream(std::istream& is)
{
	char tmp;
	while(std::isspace(is.peek()))
	{
		is.get(tmp);
	}
	return is.peek();
}

inline std::vector<std::string> file_to_lines_strip(std::string const& filename)
{
  std::ifstream inp(filename);
  if(inp.fail()) { throw std::runtime_error("Can't open file!"); }
  std::string line_buf;
  std::vector<std::string> lines;
  while(std::getline(inp, line_buf))
  {
    if(line_buf.empty()) { continue; }
    if(first_non_blank(line_buf) == '#') { continue; }
    auto begin_comment = line_buf.find('#');
    if(begin_comment != std::string::npos) { line_buf.erase(begin_comment); }
    lines.push_back(line_buf);
  }
  return lines;
}



inline constexpr double ang_to_rad(double ang) noexcept { return ang*M_PI/180.0; }
inline constexpr double rad_to_ang(double rad) noexcept { return rad*180.0/M_PI; }

inline std::size_t cell_number(std::size_t index, std::size_t orig_size) noexcept { return index / orig_size; }
inline std::size_t reduce_index(std::size_t index, std::size_t orig_size) noexcept { return index % orig_size; }

inline std::vector<std::size_t> transform_to_original_indices(std::vector<std::size_t> indices, std::size_t num_atoms) noexcept
{
  for(auto& index : indices) { index = reduce_index(index, num_atoms); }
  return indices;
}

inline double sphere_surface(double r) noexcept
{
  return 4.0 * M_PI * r * r;
}

template<typename T>
inline void append_vector(std::vector<T> const& source, std::vector<T>& dest)
{
  dest.insert(std::end(dest), std::cbegin(source), std::cend(source));
}

inline constexpr bool within_range(double a, double b, double tol)
{
  return std::abs(a-b) <= tol;
}

template<typename T>
inline bool is_subset(std::unordered_set<T> const& smaller, std::unordered_set<T> const& bigger)
{
  for(auto const& elem : smaller)
  {
    if(!bigger.contains(elem)) { return false; }
  }
  return true;
}
