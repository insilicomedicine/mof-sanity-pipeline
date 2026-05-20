// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#pragma once

#include "cif.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <vector>


namespace rdkit_wrappers
{
void determine_bonds(std::vector<atom> const& atoms, int charge=0);

std::vector<std::size_t> calc_canon_order_rdkit(std::vector<atom> const& atoms, int charge=0);

std::string xyz_stream_to_canon_smiles(std::istream *xyz_ss);
std::string atoms_to_canon_smiles(std::vector<atom> const& atoms);
}
