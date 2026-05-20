// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"
#include <optional>
#include <vector>


std::optional<std::vector<std::size_t>> match_61(std::size_t target, linker const& linker);
std::optional<std::vector<std::size_t>> match_62(std::size_t target, linker const& linker);
std::optional<std::vector<std::size_t>> match_63(std::size_t target, linker const& linker);
std::optional<std::vector<std::size_t>> match_64(std::size_t target, linker const& linker);
std::optional<std::vector<std::size_t>> match_65(std::size_t target, linker const& linker);

gene match_atom_group_3(std::size_t target, linker const& linker);

std::vector<gene> find_sites_group_1(std::vector<atom> const& atoms);
std::vector<gene> find_sites_group_2(std::vector<atom> const& atoms);
std::vector<gene> find_sites_group_3(cif_data& cif,
                                     double r_min, double r_max,
                                     double tolerance,
                                     std::string const& strategy);

bool is_group_1(int code);
bool is_group_2(int code);
bool is_group_3(int code);

void mutate_61(std::size_t target, cif_data& cif);
void mutate_62(std::size_t target, cif_data& cif);
void mutate_63(std::size_t target, cif_data& cif);
void mutate_64(std::size_t target, cif_data& cif);

void mutate_group_1(std::size_t target, int new_mutation, cif_data& cif);
void mutate_group_2(std::size_t target, int new_mutation, cif_data& cif);
void mutate_group_3(std::size_t target, int new_mutation, cif_data& cif);

