# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import pymatgen.analysis as analysis
from pymatgen.core import Structure
from pymatgen.analysis.local_env import JmolNN
from ruamel.yaml import YAML
from pymatgen.core.periodic_table import Specie

MODULE_DIR = analysis.__path__[0]


class JmolNN_update(JmolNN):
    def __init__(
        self,
        structure_to_work_with,
        tol: float = 0.45,
        min_bond_distance: float = 0.5,
    ):
        self.tol = tol
        self.min_bond_distance = min_bond_distance
        bonds_file = f"{MODULE_DIR}/bonds_jmol_ob.yaml"
        with open(bonds_file) as file:
            yaml = YAML()
            self.el_radius = yaml.load(file)

        self.structure_to_work_with = structure_to_work_with
        self.bonds = self.get_bonds()
        self.max_rad = max([bond for bonds in self.bonds.values()
                           for bond in bonds.values()]) + self.tol
        self.structure_get_neighbors = self.structure_to_work_with.get_all_neighbors(
            self.max_rad)

    def get_bonds(self):
        bonds = {specie1: {specie2: self.get_max_bond_distance(specie1.symbol, specie2.symbol)
                           for specie2 in self.structure_to_work_with.elements}
                 for specie1 in self.structure_to_work_with.elements}
        if Specie('P') and Specie('H') in self.structure_to_work_with.elements:
            bonds[Specie('P')][Specie('H')] += 1.0
            bonds[Specie('H')][Specie('P')] += 1.0
        return bonds

    def get_nn_info(self, structure: Structure, n: int):
        site = structure[n]
        min_rad = min(self.bonds[site.specie].values())
        siw = []
        for nn in self.structure_get_neighbors[n]:
            dist = nn.nn_distance
            # Confirm neighbor based on bond length specific to atom pair
            if dist <= (self.bonds[site.specie][nn.specie]) and (nn.nn_distance > self.min_bond_distance):
                weight = min_rad / dist
                siw.append(
                    {
                        "site": nn,
                        "image": self._get_image(structure, nn),
                        "weight": weight,
                        "site_index": self._get_original_site(structure, nn),
                    }
                )
        return siw