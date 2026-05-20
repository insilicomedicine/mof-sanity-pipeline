// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"
#include <array>

struct euler_angles
{
  double alpha;
  double beta;
  double gamma;
};

[[nodiscard]] std::array<double,9> calc_rotation_matrix(euler_angles const& angles);

void rotate_by_matrix(atom& a, std::array<double,9> const& rot_matrix);

void rotate_around_origin(std::vector<atom>& atoms, euler_angles const& orientation);

[[nodiscard]] euler_angles calc_euler_angles(vec3d const& orientation_from, vec3d const& normal_from,
                                             vec3d const& orientation_to, vec3d const& normal_to,
                                             [[maybe_unused]] double eps = 0.01);

void translate(std::vector<atom>& atoms, vec3d const& translation);

std::vector<atom> kabsch_algorithm(std::vector<atom> const& target, std::vector<std::size_t> const& target_indices,
                                   std::vector<atom> const& ref, std::vector<std::size_t> const& ref_indices);
