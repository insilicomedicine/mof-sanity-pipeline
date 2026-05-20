// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "utils.hpp"
#include "rdkit_wrappers.hpp"

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/DetermineBonds/DetermineBonds.h>

#include <sstream>
#include <memory>

namespace rdkit_wrappers
{
std::vector<unsigned int> CanonicalRankAtoms(RDKit::ROMol const& mol,
                                             bool breakTies = true, bool includeChirality = true,
                                             bool includeIsotopes = true, bool includeAtomMaps = true,
                                             bool includeChiralPresence = false)
{
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  const bool includeStereoGroups = true;

  RDKit::Canon::rankMolAtoms(mol, ranks, breakTies,
                             includeChirality, includeIsotopes,
                             includeAtomMaps, includeChiralPresence,
                             includeStereoGroups);
  return ranks;
}

void determine_bonds(RDKit::RWMol& mol, int charge)
{
  bool const useHueckel = false;
  double const covFactor = 1.3;
  bool const useVdw = true;
  bool allowChargedFragments = true;
  bool embedChiral = true;
  bool useAtomMap = false;
  RDKit::determineBonds(mol, useHueckel, charge, covFactor,
                        allowChargedFragments, embedChiral,
                        useAtomMap, useVdw);
}

void determine_bonds(std::vector<atom> const& atoms, int charge)
{
  std::stringstream ss;
  print_xyz(atoms, &ss);

  std::unique_ptr<RDKit::RWMol> mol_uptr(RDKit::XYZDataStreamToMol(ss));
  determine_bonds(*mol_uptr, charge);
}

void determine_connectivity(RDKit::RWMol& mol)
{
  bool const useHueckel = false;
  double const covFactor = 1.3;
  bool const useVdw = true;
  int const charge = 0;
  RDKit::determineConnectivity(mol, useHueckel, charge, covFactor, useVdw);
}

std::vector<std::size_t> ranks_to_order(std::vector<unsigned int> const& ranks)
{
  std::map<unsigned int, std::size_t> ranks_with_indices;
  for(std::size_t i=0; i<ranks.size(); ++i) { ranks_with_indices[ranks[i]] = i; }
  std::vector<std::size_t> canon_order;
  for(auto const& [rank,index] : ranks_with_indices)
  {
    canon_order.push_back(index);
  }
  return canon_order;
}

std::vector<std::size_t> tmp_func(std::vector<unsigned int> const& ranks)
{
  std::vector<std::size_t> res;
  for(auto rank : ranks)
  {
    res.push_back(static_cast<std::size_t>(rank));
  }
  return res;
}

std::vector<std::size_t> calc_canon_order_rdkit(std::vector<atom> const& atoms, [[maybe_unused]] int charge)
{
  std::stringstream ss;
  print_xyz(atoms, &ss);

  std::unique_ptr<RDKit::RWMol> mol_uptr(RDKit::XYZDataStreamToMol(ss));

  determine_connectivity(*mol_uptr);

  auto canon_ranks = CanonicalRankAtoms(*mol_uptr);
  auto canon_order = tmp_func(canon_ranks);
  return canon_order;
}

std::string xyz_stream_to_canon_smiles(std::istream *xyz_ss)
{
  std::unique_ptr<RDKit::ROMol> mol_uptr(RDKit::XYZDataStreamToMol(*xyz_ss));
  return MolToSmiles(*mol_uptr);
}

std::string atoms_to_canon_smiles(std::vector<atom> const& atoms)
{
  std::stringstream ss;
  print_xyz(atoms, &ss);
  return xyz_stream_to_canon_smiles(&ss);
}

}
