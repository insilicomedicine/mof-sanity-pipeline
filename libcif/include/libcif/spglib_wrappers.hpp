// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"

namespace spglib_wrappers
{

void primitivize_cell(cif_data& cif, double tolerance = 0.25);

}
