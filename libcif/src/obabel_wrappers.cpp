// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Saudi Aramco -- MIT License.
// See repository LICENSE for full terms.
#include "utils.hpp"
#include "obabel_wrappers.hpp"
#include <openbabel3/openbabel/obconversion.h>

namespace obabel_wrappers
{

std::string xyz_stream_to_canon_smiles(std::istream* xyz_ss)
{
  std::ostringstream os;
  OpenBabel::OBConversion conv(xyz_ss, &os);
  conv.SetInAndOutFormats("xyz", "can");
  conv.Convert();
  return os.str();
}

std::string atoms_to_canon_smiles(std::vector<atom> const& atoms)
{
  std::stringstream ss;
  print_xyz(atoms, &ss);
  auto smiles = xyz_stream_to_canon_smiles(&ss);
  return smiles;
}

}
