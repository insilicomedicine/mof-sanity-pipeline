// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"
#include <vector>


namespace obabel_wrappers
{

std::string xyz_stream_to_canon_smiles(std::istream* xyz_ss);
std::string atoms_to_canon_smiles(std::vector<atom> const& atoms);

}
