// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "transform_ops.hpp"
#include "cif_ops.hpp"
#include "linalg.hpp"

[[nodiscard]] std::array<double,9> calc_rotation_matrix(euler_angles const& angles)
{
  double const cosa = std::cos(angles.alpha);
  double const cosb = std::cos(angles.beta);
  double const cosg = std::cos(angles.gamma);
  double const sina = std::sin(angles.alpha);
  double const sinb = std::sin(angles.beta);
  double const sing = std::sin(angles.gamma);

  std::array<double,9> rot_matrix;
  rot_matrix[0] = cosa*cosg - sina*cosb*sing;
  rot_matrix[1] = -cosa*sing - sina*cosb*cosg;
  rot_matrix[2] = sina*sinb;
  rot_matrix[3] = sina*cosg + cosa*cosb*sing;
  rot_matrix[4] = -sina*sing + cosa*cosb*cosg;
  rot_matrix[5] = -cosa*sinb;
  rot_matrix[6] = sinb*sing;
  rot_matrix[7] = sinb*cosg;
  rot_matrix[8] = cosb;

  return rot_matrix;
}

void rotate_by_matrix(atom& a, std::array<double,9> const& rot_matrix)
{
  std::array<double,3> new_coord{0,0,0};
  for(int k=0; k<3; ++k)
  {
    for(int l=0; l<3; ++l)
    {
      new_coord[k] += a.coord[l] * rot_matrix[k*3+l];
    }
  }
  for(int k=0; k<3; ++k) { a.coord[k] = new_coord[k]; }
}

void rotate_around_origin(std::vector<atom>& atoms, euler_angles const& orientation)
{
  if(orientation.alpha == 0 and orientation.beta == 0 and orientation.gamma == 0) { std::cerr << "NOT ROTATING\n"; return; }
  auto rot_matrix = calc_rotation_matrix(orientation);
  for(auto& atom : atoms)
  {
    rotate_by_matrix(atom, rot_matrix);
  }
}

[[nodiscard]] euler_angles calc_euler_angles(vec3d const& orientation_from, vec3d const& normal_from,
                                             vec3d const& orientation_to, vec3d const& normal_to,
                                             [[maybe_unused]] double eps)
{
  auto aux_from_x = orientation_from;
  auto aux_from_y = vector_product(aux_from_x, normal_from);
  auto aux_from_z = vector_product(orientation_from, aux_from_y);

  auto aux_to_x = orientation_to;
  auto aux_to_y = vector_product(aux_to_x, normal_to);
  auto aux_to_z = vector_product(orientation_to, aux_to_y);

  euler_angles angles;

  double const m33 = aux_from_x[2]*aux_to_x[2] +
                     aux_from_y[2]*aux_to_y[2] +
                     aux_from_z[2]*aux_to_z[2];


  if(std::abs(m33) < 1.0 - 1e-10)
  {
    double const m31 = aux_from_x[0]*aux_to_x[2] + 
      aux_from_y[0] * aux_to_y[2] + 
      aux_from_z[0] * aux_to_z[2];
    double const m32 = aux_from_x[1] * aux_to_x[2] +
      aux_from_y[1] * aux_to_y[2] +
      aux_from_z[1] * aux_to_z[2];
    double const m13 = aux_from_x[2] * aux_to_x[0] +
      aux_from_y[2] * aux_to_y[0] +
      aux_from_z[2] * aux_to_z[0];
    double const m23 = aux_from_x[2] * aux_to_x[1] +
      aux_from_y[2] * aux_to_y[1] +
      aux_from_z[2] * aux_to_z[1];

    angles.beta = std::acos(m33);
    angles.gamma = std::atan2(m31, m32);
    angles.alpha = std::atan2(m13, -m23);
  }
  else
  {
    double const m11 = aux_from_x[0] * aux_to_x[0] +
      aux_from_y[0] * aux_to_y[0] +
      aux_from_z[0] * aux_to_z[0];
    double const m21 = aux_from_x[0] * aux_to_x[1] +
      aux_from_y[0] * aux_to_y[1] +
      aux_from_z[0] * aux_to_z[1];

    angles.alpha = std::atan2(m21, m11);
    angles.beta = m33<0 ? M_PI : 0;
    angles.gamma = 0;
  }
  
  return angles;
}

void translate(std::vector<atom>& atoms, vec3d const& translation)
{
  for(auto& atom : atoms) { atom.coord += translation; }
}

vec3d calc_centroid(std::vector<atom> const& atoms)
{
  vec3d center {0,0,0};
  for(auto const& atom : atoms) { center += atom.coord; }
  center /= atoms.size();
  return center;
}

std::vector<atom> kabsch_algorithm(std::vector<atom> const& target, std::vector<std::size_t> const& target_indices,
                                   std::vector<atom> const& ref, std::vector<std::size_t> const& ref_indices)
{
  if(target_indices.size() != ref_indices.size()) { throw std::runtime_error("Anchor indices have different size!"); }

  auto target_anchors = indices_to_atoms(target_indices, target);
  auto ref_anchors = indices_to_atoms(ref_indices, ref);

  auto target_centroid = calc_centroid(target_anchors);
  auto ref_centroid = calc_centroid(ref_anchors);

  translate(target_anchors, -target_centroid);
  translate(ref_anchors, -ref_centroid);

  matrix33 cov_matrix {};
  for(std::size_t i=0; i<target_anchors.size(); ++i)
  {
    cov_matrix += outer_product(ref_anchors[i].coord, target_anchors[i].coord);
  }

  auto rot_matrix = find_rotation_svd(cov_matrix);

  auto translation = ref_centroid - multiply(rot_matrix, target_centroid);

  std::vector<atom> res = target;
  for(auto& atom : res)
  {
    atom.coord = multiply(rot_matrix, atom.coord) + translation;
  }

  return res;
}
