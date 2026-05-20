// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

template<std::size_t N>
class vec
{
private:
  std::array<double,N> data {};

public:
  vec() = default;
  vec(double val1, double val2, double val3) : data{val1, val2, val3} {}
  vec(const vec& other) : data(other.data) {}
  vec& operator=(const vec& other)
  {
    this->data = other.data;
    return *this;
  }
  ~vec() = default;

  explicit vec(std::array<double,N> arr) : data(std::move(arr)) {}

  double& operator[](std::size_t index)
  {
    return data[index];
  }
  double const& operator[](std::size_t index) const
  {
    return data[index];
  }

  vec& operator+=(double val)
  {
    for(std::size_t i=0; i<N; ++i)
    {
      this->data[i] += val;
    }
    return *this;
  }
  vec& operator-=(double val)
  {
    for(std::size_t i=0; i<N; ++i)
    {
      this->data[i] -= val;
    }
    return *this;
  }
  vec& operator*=(double val)
  {
    for(std::size_t i=0; i<N; ++i)
    {
      this->data[i] *= val;
    }
    return *this;
  }
  vec& operator/=(double val)
  {
    for(std::size_t i=0; i<N; ++i)
    {
      this->data[i] /= val;
    }
    return *this;
  }

  vec& operator+=(vec<N> const& other)
  {
    for(std::size_t k=0; k<N; ++k)
    {
      this->data[k] += other.data[k];
    }
    return *this;
  }
  vec& operator-=(vec<N> const& other)
  {
    for(std::size_t k=0; k<N; ++k)
    {
      this->data[k] -= other.data[k];
    }
    return *this;
  }

  vec operator+() const
  {
    return *this;
  }
  vec operator-() const
  {
    vec<N> res;
    for(std::size_t i=0; i<N; ++i)
    {
      res.data[i] = -data[i];
    }
    return res;
  }

  std::array<double,N> to_array()
  {
    return this->data;
  }
};

template<std::size_t N, std::size_t M>
class matrix
{
private:
  static constexpr std::size_t flat_size = N*M;
  std::array<double,flat_size> data {}; // row-major

public:
  double& operator()(std::size_t row, std::size_t col)
  {
    return data[row*M + col];
  }
  double const& operator()(std::size_t row, std::size_t col) const
  {
    return data[row*M + col];
  }

  matrix& operator *=(double val)
  {
    for(std::size_t i=0; i<flat_size; ++i) { this->data[i] *= val; }
    return *this;
  }
  matrix& operator /=(double val)
  {
    for(std::size_t i=0; i<flat_size; ++i) { this->data[i] /= val; }
    return *this;
  }

  matrix& operator +=(matrix<N,M> const& other)
  {
    for(std::size_t i=0; i<flat_size; ++i) { this->data[i] += other.data[i]; }
    return *this;
  }
  matrix& operator -=(matrix<N,M> const& other)
  {
    for(std::size_t i=0; i<flat_size; ++i) { this->data[i] -= other.data[i]; }
    return *this;
  }

  matrix operator+() const { return *this; }
  matrix operator-() const
  {
    matrix<N,M> result;
    for(std::size_t i=0; i<flat_size; ++i) { result[i] = -data[i]; }
    return result;
  }
};

using vec3d = vec<3>;
using matrix33 = matrix<3,3>;

template<std::size_t N>
vec<N> operator+(vec<N> v, double val)
{
  v += val;
  return v;
}

template<std::size_t N>
vec<N> operator-(vec<N> v, double val)
{
  v -= val;
  return v;
}

template<std::size_t N>
vec<N> operator*(vec<N> v, double val)
{
  v *= val;
  return v;
}

template<std::size_t N>
vec<N> operator*(double val, vec<N> v)
{
  v *= val;
  return v;
}

template<std::size_t N>
vec<N> operator/(vec<N> v, double val)
{
  v /= val;
  return v;
}

template<std::size_t N>
vec<N> operator+(vec<N> lhs, vec<N> const& rhs)
{
  lhs += rhs;
  return lhs;
}

template<std::size_t N>
vec<N> operator-(vec<N> lhs, vec<N> const& rhs)
{
  lhs -= rhs;
  return lhs;
}

template<std::size_t N>
inline double scalar_product(vec<N> const& a, vec<N> const& b)
{
  double res = 0;
  for(std::size_t k=0; k<N; ++k)
  {
    res += a[k] * b[k];
  }
  return res;
}

template<std::size_t N>
[[nodiscard]] inline double norm(vec<N> const& v)
{
  return std::sqrt(scalar_product(v, v));
}

template<std::size_t N>
inline void normalize(vec<N>& v)
{
  v /= norm(v);
}

template<std::size_t N>
inline void normalize(vec<N>& v, double norm)
{
  v /= norm;
}

template<std::size_t N>
inline double calc_norm_and_normalize(vec<N>& v)
{
  double n = norm(v);
  if(n == 0) { std::cerr << "WARNING!! CAN'T NORMALIZE ZERO VECTOR\n"; return 0; }
  v /= n;
  return n;
}

[[nodiscard]] inline vec3d vector_product(vec3d const& a, vec3d const& b)
{
  return
  {
     a[1]*b[2] - a[2]*b[1],
    -a[0]*b[2] + a[2]*b[0],
     a[0]*b[1] - a[1]*b[0]
  };
}

inline vec3d rotate_around_axis(vec3d const& vec, vec3d const& axis, double angle)
{
  // Using Rodrigues formula
  // v_rot = v*cos(theta) + (k x v)*sin(theta) + k*(kv)*(1-cos(theta))
  // where v -- original, k -- unit axis vector

  // First term
  vec3d rotated_vec = vec*std::cos(angle);

  // Second term
  auto second_term = vector_product(vec, axis) * std::sin(angle);

  // Third term
  auto third_term = axis * scalar_product(vec, axis) * (1 - std::cos(angle));

  rotated_vec += second_term;
  rotated_vec += third_term;

  return rotated_vec;
}

template<std::size_t N>
inline double angle(vec<N> const& vec1, vec<N> const& vec2)
{
  return std::acos(scalar_product(vec1, vec2) / norm(vec1) / norm(vec2));
}

template<std::size_t N, std::size_t M>
inline matrix<N,M> outer_product(vec<N> const& rhs, vec<M> const& lhs)
{
  matrix<N,M> res;
  for(std::size_t i=0; i<N; ++i)
  {
    for(std::size_t j=0; j<M; ++j)
    {
      res(i,j) = rhs[i] * lhs[j];
    }
  }
  return res;
}



// ***********************************
// MATRIX OPERATIONS
// ***********************************

template<std::size_t N>
matrix<N,N> make_unit_matrix()
{
  matrix<N,N> unit;
  for(std::size_t i=0; i<N; ++i) { unit(i,i) = 1; }
  return unit;
}

template<std::size_t N, std::size_t M>
matrix<N,M> operator*(matrix<N,M> m, double val)
{
  m *= val;
  return m;
}

template<std::size_t N, std::size_t M>
matrix<N,M> operator*(double val, matrix<N,M> m)
{
  m *= val;
  return m;
}

template<std::size_t N, std::size_t M>
matrix<N,M> operator+(matrix<N,M> rhs, matrix<N,M> const& lhs)
{
  rhs += lhs;
  return rhs;
}

template<std::size_t N, std::size_t M>
matrix<N,M> operator-(matrix<N,M> rhs, matrix<N,M> const& lhs)
{
  rhs -= lhs;
  return rhs;
}

template<std::size_t N, std::size_t M>
inline matrix<M,N> transpose(matrix<N, M> const& m)
{
  matrix<M,N> res;
  for(std::size_t i=0; i<N; ++i)
  {
    for(std::size_t j=0; j<M; ++j)
    {
      res(j,i) = m(i,j);
    }
  }
  return res;
}

inline double determinant(matrix33 const& m)
{
  return m(0,0) * (m(1,1)*m(2,2) - m(1,2)*m(2,1)) -
         m(0,1) * (m(1,0)*m(2,2) - m(1,2)*m(2,0)) +
         m(0,2) * (m(1,0)*m(2,1) - m(1,1)*m(2,0));
}

inline matrix33 inverse(matrix33 const& m)
{
  matrix33 inv;
  inv(0,0) = m(1,1)*m(2,2) - m(1,2)*m(2,1);
  inv(0,1) = m(0,2)*m(2,1) - m(0,1)*m(2,2);
  inv(0,2) = m(0,1)*m(1,2) - m(0,2)*m(1,1);
  
  inv(1,0) = m(1,2)*m(2,0) - m(1,0)*m(2,2);
  inv(1,1) = m(0,0)*m(2,2) - m(0,2)*m(2,0);
  inv(1,2) = m(0,2)*m(1,0) - m(0,0)*m(1,2);
  
  inv(2,0) = m(1,0)*m(2,1) - m(1,1)*m(2,0);
  inv(2,1) = m(0,1)*m(2,0) - m(0,0)*m(2,1);
  inv(2,2) = m(0,0)*m(1,1) - m(0,1)*m(1,0);

  double det = determinant(m);
  if(std::abs(det) < 1e-12) { return make_unit_matrix<3>(); }
  inv /= det;

  return inv;
}

template<std::size_t N1, std::size_t M1, std::size_t N2, std::size_t M2>
inline matrix<N1,M2> multiply(matrix<N1,M1> const& rhs, matrix<N2,M2> const& lhs)
{
  static_assert(M1==N2, "wrong matrix size!");

  matrix<N1,M2> res;
  for(std::size_t i=0; i<N1; ++i)
  {
    for(std::size_t j=0; j<M2; ++j)
    {
      double sum = 0;
      for(std::size_t k=0; k<M1; ++k)
      {
        sum += rhs(i,k) * lhs(k,j);
      }
      res(i,j) = sum;
    }
  }
  return res;
}

template<std::size_t N, std::size_t M>
inline vec<N> multiply(matrix<N,M> const& m, vec<M> const& v)
{
  vec<N> res;
  for(std::size_t i=0; i<N; ++i)
  {
    double sum = 0;
    for(std::size_t k=0; k<M; ++k)
    {
      sum += m(i,k) * v[k];
    }
    res[i] = sum;
  }
  return res;
}

inline matrix33 find_rotation_polar(matrix33 const& m, int max_it = 10, double tolerance = 1e-9)
{
  matrix33 r = m;

  for(int n=0; n<max_it; ++n)
  {
    auto r_trans = transpose(r);
    auto r_trans_inv = inverse(r_trans);
    
    auto r_new = 0.5 * (r + r_trans_inv);

    // Optional check for early convergence
    matrix33 ortho_check = multiply(transpose(r_new), r_new);
    matrix33 unit = make_unit_matrix<3>();

    double max_diff = 0;
    for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
      {
        max_diff = std::max(max_diff, std::abs(ortho_check(i,j) - unit(i,j)));
      }
    }

    r = r_new;

    if(max_diff<tolerance && n>2) { break; }
  }

  // Proper rotation
  double det = determinant(r);
  // if det<0 flip third column
  if(det<0)
  {
    for(int i=0; i<3; ++i)
    {
      r(i,2) *= -1;
    }
  }

  return r;
}

// ---------------------------------------------------------
// JACOBI SVD IMPLEMENTATION FOR 3x3
// ---------------------------------------------------------

// Helper to perform a Jacobi rotation on a 3x3 matrix
inline void jacobi_rotate(matrix33& A, matrix33& R, int p, int q)
{
  if (std::abs(A(p,q)) < 1e-12) return;

  double d = (A(q,q) - A(p,p)) / (2.0 * A(p,q));
  double t = 1.0 / (std::abs(d) + std::sqrt(d*d + 1.0));
  if (d < 0) t = -t;

  double c = 1.0 / std::sqrt(t*t + 1.0);
  double s = t * c;

  // Update A elements manually for symmetry preservation
  double App = A(p,p);
  double Aqq = A(q,q);
  double Apq = A(p,q);

  A(p,p) = c*c*App - 2.0*s*c*Apq + s*s*Aqq;
  A(q,q) = s*s*App + 2.0*s*c*Apq + c*c*Aqq;
  A(p,q) = 0.0;
  A(q,p) = 0.0;

  for (int i = 0; i < 3; ++i) {
    if (i != p && i != q) {
      double Api = A(p,i);
      double Aqi = A(q,i);
      A(p,i) = c*Api - s*Aqi;
      A(q,i) = s*Api + c*Aqi;
      A(i,p) = A(p,i);
      A(i,q) = A(q,i);
    }
  }

  // Update V matrix
  for (int i = 0; i < 3; ++i) {
    double Rip = R(i,p);
    double Riq = R(i,q);
    R(i,p) = c*Rip - s*Riq;
    R(i,q) = s*Rip + c*Riq;
  }
}

// Complete SVD-based rotation finder
inline matrix33 find_rotation_svd(matrix33 const& H)
{
  // 1. Construct S = H^T * H
  matrix33 H_t = transpose(H);
  matrix33 S = multiply(H_t, H);

  // 2. Diagonalize S to find V
  matrix33 V = make_unit_matrix<3>();
  
  for (int iter = 0; iter < 30; ++iter) {
    int p = 0, q = 1;
    double max_val = std::abs(S(0,1));
    if(std::abs(S(0,2)) > max_val) { p=0; q=2; max_val = std::abs(S(0,2)); }
    if(std::abs(S(1,2)) > max_val) { p=1; q=2; max_val = std::abs(S(1,2)); }

    if (max_val < 1e-12) break; 
    jacobi_rotate(S, V, p, q);
  }

  // 3. Compute U = H * V
  matrix33 U = multiply(H, V);

  // 4. Normalize columns of U and handle Singularity (Flat Case)
  int singular_col = -1; 
  
  for (int col = 0; col < 3; ++col) {
    double mag = 0.0;
    for (int row = 0; row < 3; ++row) {
      mag += U(row, col) * U(row, col);
    }
    mag = std::sqrt(mag);

    if (mag > 1e-12) {
      for (int row = 0; row < 3; ++row) {
        U(row, col) /= mag;
      }
    } else {
       singular_col = col; // Found a zero column (flat dimension)
    }
  }

  // RECONSTRUCTION LOGIC:
  // If we have a singular column, we construct it as the cross product 
  // of the other two valid columns to ensure orthonormality.
  if (singular_col != -1) {
      int c1 = (singular_col + 1) % 3;
      int c2 = (singular_col + 2) % 3;
      
      vec3d v1 { U(0, c1), U(1, c1), U(2, c1) };
      vec3d v2 { U(0, c2), U(1, c2), U(2, c2) };
      
      vec3d v_sing = vector_product(v1, v2);
      
      // If v_sing is also zero (linear molecule), just pick arbitrary axis
      if(norm(v_sing) < 1e-12) {
          // Fallback for linearity (rare for Kabsch but possible)
          v_sing = {1,0,0}; 
          // If v1 was {1,0,0}, pick {0,1,0}
          if(std::abs(v1[0]) > 0.9) v_sing = {0,1,0};
      } else {
          normalize(v_sing);
      }

      U(0, singular_col) = v_sing[0];
      U(1, singular_col) = v_sing[1];
      U(2, singular_col) = v_sing[2];
  }

  // 5. Calculate R = U * V^T
  matrix33 R = multiply(U, transpose(V));

  // 6. Check determinant to ensure a proper rotation (determinant +1)
  // If determinant is -1, it's a reflection. We must flip the axis corresponding
  // to the smallest singular value (which is likely the one we just reconstructed).
  double det = determinant(R);
  
  if (det < 0) {
    // Find column in S corresponding to smallest eigenvalue
    int min_idx = 0;
    double min_val = S(0,0);
    if(S(1,1) < min_val) { min_val = S(1,1); min_idx = 1; }
    if(S(2,2) < min_val) { min_val = S(2,2); min_idx = 2; }
    
    // Flip that column in V (effectively U * Flip * V^T)
    // We can just construct a correction matrix
    matrix33 Correction = make_unit_matrix<3>();
    Correction(min_idx, min_idx) = -1.0;
    
    R = multiply(multiply(U, Correction), transpose(V));
  }
  
  return R;
}
