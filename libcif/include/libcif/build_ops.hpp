// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "building_blocks.hpp"
#include "cif.hpp"

[[nodiscard]] cif_data build_hum(building_blocks::topology const& topology);
[[nodiscard]] cif_data build_999(building_blocks::topology const& topology, double const c = 5.0);

[[nodiscard]] building_blocks::node_type read_xyz_to_node(std::string const& filename);
[[nodiscard]] building_blocks::edge_type read_xyz_to_edge(std::string const& filename);
