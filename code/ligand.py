# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import numpy as np
from openbabel import openbabel
from rdkit import Chem
from dataclasses import dataclass
from rdkit.Chem import rdDetermineBonds
import networkx as nx


from openbabel import pybel
import warnings

warnings.simplefilter('ignore')

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
pybel.ob.obErrorLog.StopLogging()


import logging
logger = logging.getLogger(__name__)

SMARTS_ANION_DICT = {
    "[F,Cl,Br,I;X0]" : ("Halogenide anion [Hal]-" , -1, [0]),
    "[O,S,Se,Te;X0]" : ('Oxygen/S/Se/Te Lone Atom 2-', -2, [0]),
    "[CX1]~[NX1]" : ("Cyanide or isocyanide [CN]-", -1, (0,1)),
    "[OX1]~[CX2]~[NX1]" : ("Isocyanate [OCN]-", -1, (0,2)),
    "[CX1]~[NX2]~[OX1]" : ("Cyanate [CNO]-", -1, (0,2)),
    "[SX1]~[CX2]~[NX1]" : ("Isothiocyanate [SCN]-", -1, (0,2)),
    "[Se;X1]~[CX2]~[NX1]" : ("Isoselenocyanate [SeCN]-", -1, (0,2)),
    "[CX1]~[NX2]~[SX1]" : ("Thiocyanate [CNS]-", -1, (0,2)),
    "[CX1]~[NX2]~[Se;X1]" : ("Selenocyanate [CNSe]-", -1, (0,2)),
    "[OX1]~[#1;X1]" : ("Hydroxyl [OH]-", -1, [0]),
    "[OX1]~[PX4](~[OX1])(~[OX1])~[OX1]" : ("Phosphate anion [PO4]3-", -3, (0,2,3,4) ),
    "[#1;X1]~[PX4](~[OX1])(~[OX1])~[OX1]" : ("Phosphite anion [HPO3]2-", -2, (2,3,4) ),
    "[OX1]~[PX4](~[OX1])(~[OX1])~[OX2]~[PX4](~[OX1])(~[OX1])~[OX1]" : ("Diphosphate anion [P2O7]4-", -4, (0,2,3,6,7,8) ),
    "[OX1]~[SX4](~[OX1])(~[OX1])~[OX1]" : ("Sulfate anion [SO4]2-", -2, (0,2,3,4) ),
    "[OX1]~[SX4](~[SX1])(~[OX1])~[OX1]" : ("Thiosulfate anion [S=SO3]2-", -2, (0,3,4) ),
    "[OX1]~[SX3](~[OX1])~[OX1]" : ("Sulfite anion [SO3]2-", -2, (0,2,3) ),
    "[OX1]~[NX3](~[OX1])~[OX1]" : ("Nitrate anion [NO3]-", -1, (0,2,3) ),
    "[OX1]~[NX2]~[OX1]" : ("Nitrite anion [NO2]-", -1, (0,2) ),
    "[NX1]~[NX2]~[NX1]" : ("Azide anion [N3]-", -1, (0,2) ),
    "[NX1]~[NX1]" : ("-N#N- [N2]-", -1, (0,1) ),
    "[#1;X1]~[#8;X2]~[#1;X1]" : ("Water", 0, ()),
    "[OX1]~[CX2]~[OX1]" : ("CO2", 0, ()),
    "[C;X1]~[O;X1]" : ("CO", 0, ()),
    "[N;X1]~[N;X1]" : ("N2", 0, ()),
    "[#1;X1]~[#1;X1]" : ("H2", 0, ()),
    "[#8;X1]~[#8;X1]" : ("O2", 0, ()),
    "[#1;X1]~[N;X3](~[#1;X1])~[#1;X1]" : ("Ammonia", 0, ()),
    "[O;X1]~[Cl,Br,I;X4](~[O;X1])(~[O;X1])~[O;X1]" : ("Perhalogenate anion [HalO4]-", -1, (0,2,3,4) ),
    "[O;X1]~[Cl,Br,I;X3](~[O;X1])~[O;X1]" : ("Halogenate anion [HalO3]-", -1, (0,2,3) ),
    "[O;X1]~[Cl,Br,I;X2]~[O;X1]" : ("Halogenite anion [HalO2]-", -1, (0,2) ),
    "[O;X1]~[Cl,Br,I;X1]" : ("Hypohalogenite anion [HalO]-", -1, [0] ),
    "[Cl,Br,I;X1]~[Br,I;X2]~[Cl,Br,I;X1]" : ("[Hal_X Hal_Y_2] -", -1, [0]),
    "[Cl,Br,I;X1]~[Cl,Br,I;X1]" : ("[Hal_2]", 0, () ),
    "[#1;X1]~[Cl,Br,I;X1]" : ("[H Hal]", 0, () )
}

class Anion:
    def __init__(self, smarts, name, total_charge, charged_atoms):
        self.name = name
        self.smarts = smarts
        self.mol = Chem.MolFromSmarts(smarts)
        self.total_charge = total_charge
        self.charged_atoms = charged_atoms
        self.N = self.mol.GetNumAtoms()

Anions = [ Anion(smarts, *data) for smarts, data in SMARTS_ANION_DICT.items() ]

@dataclass
class ValenceAtom:
    symbol : str
    standart_valence : int
    may_be_higher : bool = True
    may_be_lower : bool = False
    multiplier : int = 1

    def get_adding_charge(self, valence):
        if ( valence>self.standart_valence and not self.may_be_higher ) or \
                (valence<self.standart_valence and not self.may_be_lower):
                raise Exception(f"ERROR {self.symbol} {self.standart_valence} {valence}")
        return  (valence-self.standart_valence)*self.multiplier


ValenceAtoms = {
    'C'  : ValenceAtom('C', 4, False, True, -1),
    'N'  : ValenceAtom('N', 3, True, True),
    'O'  : ValenceAtom('O', 2, False, True),
    'F'  : ValenceAtom('F', 1, False, True),
    'Cl'  : ValenceAtom('Cl', 1, False, True),
    'Br'  : ValenceAtom('Br', 1, False, True),
    'I'  : ValenceAtom('I', 1, False, True),
    'H'  : ValenceAtom('H', 1, False, True, -1),
    'Si' : ValenceAtom('Si', 4, True, False, -1),
    'Ge' : ValenceAtom('Ge', 4, True, False, -1),
    'S'  : ValenceAtom('S', 2, False, True),
    'P'  : ValenceAtom('P', 5, True, True)
}

class GroupedCharge:
    def __init__(self, idxs, charge):
        self.idxs = idxs
        self.charge = charge
        self.metal_counter = 0

    def detect_metals_charges_and_free_charge(self, g, metal_idxs, metal_charges):
        n = 0
        metal_n_counter = {}
        for metal_idx in metal_idxs:
            metal_n_counter[metal_idx] = len(set([s.index for s in g.get_connected_sites(metal_idx) if s.jimage==(0,0,0)]).intersection(set(self.idxs)))
            n += metal_n_counter[metal_idx]
        if n==0:
            return self.charge
        for metal_idx in metal_idxs:
            metal_charges[metal_idx] += -1*metal_n_counter[metal_idx]*self.charge/n
        return 0

class Ligand:

    @staticmethod
    def detect_ligand(idx, g, visited_sites, ligand, current_jimage, idx_in_ligands):
        if np.sum(np.array(current_jimage))>6:
            raise Exception("Infinite ligand detected")
        if [idx, *current_jimage] in visited_sites:
            return
        if idx in idx_in_ligands:
            return
        ligand.append([idx, *current_jimage])
        visited_sites.append([idx, *current_jimage])
        for site in g.get_connected_sites(idx):
            if not site.site.specie.is_metal:
                jimage = tuple(np.array(current_jimage) + np.array(site.jimage))
                Ligand.detect_ligand(site.index, g, visited_sites, ligand, jimage, idx_in_ligands)

    def gen_xyz_string(self):
        N = len(self.symbols)
        xyzstr = f"{N}\n\n"
        for symbol, coord in zip(self.symbols, self.coords):
            xyzstr += f"{symbol:10} {coord[0]:10.6f} {coord[1]:10.6f} {coord[2]:10.6f}\n"
        return xyzstr

    def get_mol_string(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats('xyz', "sdf")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, self.gen_xyz_string())
        outmol = obConversion.WriteString(mol)
        return outmol

    def gen_rdkit(self):
        molblock = self.get_mol_string()
        self.mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
        self.graph = nx.Graph()
        for a in self.mol.GetAtoms():
            self.graph.add_node(a.GetIdx(), symbol=a.GetSymbol())
        for b in self.mol.GetBonds():
            self.graph.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bondtype=b.GetBondType())
        self.graph_hashes = {"molecule" : nx.weisfeiler_lehman_graph_hash(self.graph, node_attr='symbol'),
                             "undecorated" : nx.weisfeiler_lehman_graph_hash(self.graph),
                             "with_bond_info" : nx.weisfeiler_lehman_graph_hash(self.graph,
                                                                                node_attr='symbol',
                                                                                edge_attr='bondtype')}

    def get_total_charge(self, marked_as_visited):
        charge = 0
        for idx, atom in enumerate(self.mol.GetAtoms()):
            if idx not in marked_as_visited:
                charge +=  ValenceAtoms[atom.GetSymbol()].get_adding_charge(atom.GetTotalValence())
        self.charge = charge
        return charge

    def sanitize(self):
        charge = 0
        marked_as_visited = []
        match_para_cRcR = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6;X3;v3;R1;r6]-[*;R1;r6]=[*;R1;r6]-[#6;X3;v3;R1;r6]'))
        for match in match_para_cRcR:
            self.mol.GetBondBetweenAtoms(match[0],match[1]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[2],match[3]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.SINGLE)
            marked_as_visited += list(match)
        matches_cRcR = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6;X3;H1;D3;v3]~[#6;X3;H1;D3;v3]'))
        for match in matches_cRcR:
            self.mol.GetBondBetweenAtoms(match[0],match[1]).SetBondType(Chem.BondType.DOUBLE)
            marked_as_visited += list(match)
        matches_bad_c_in_ring = self.mol.GetSubstructMatches(Chem.MolFromSmarts("[#6;r5;R1;X3;v3]-[*;r5;R1]=[#7;r5;R1;X2;v3]"))
        first_idx = []
        for match in matches_bad_c_in_ring:
            if match[0] not in first_idx:
                first_idx.append(match[0])
                self.mol.GetBondBetweenAtoms(match[0],match[1]).SetBondType(Chem.BondType.DOUBLE)
                self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.SINGLE)
                self.mol.GetAtomWithIdx(match[2]).SetFormalCharge(-1)
                marked_as_visited += list(match)
                charge -= 1

        matches_bad_other_c_in_ring = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#7&X2&v2,#6&X3&v3;R1;r5]-[#7&X2&v3,#6&X3&v4;R1;r5]=[#7&X2&v3,#6&34&v4;R1;r5]-[#7&X2&v2,#6&X3&v3;R1;r5]'))
        for match in matches_bad_other_c_in_ring:
            self.mol.GetBondBetweenAtoms(match[0],match[1]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetBondBetweenAtoms(match[2],match[3]).SetBondType(Chem.BondType.DOUBLE)
            marked_as_visited += list(match)
        matches_bad_n_in_ring = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#7;R1;r5;X2;v2]'))
        for match in matches_bad_n_in_ring:
            if match[0] not in marked_as_visited:
                self.mol.GetAtomWithIdx(match[0]).SetFormalCharge(-1)
                marked_as_visited += list(match)
                charge -= 1

        match_SO3 = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#8;X1]~[#16;X4](~[#8;X1])~[#8;X1]'))
        for match in match_SO3:
            marked_as_visited += list(match)
            charge -= 1

        match_NO2 = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]~[#7;X3](~[#8;X1])~[#8;X1]'))
        for match in match_NO2:
            marked_as_visited += list(match[1:])
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetAtomWithIdx(match[1]).SetFormalCharge(+1)
            self.mol.GetAtomWithIdx(match[3]).SetFormalCharge(-1)

        match_NO3m = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#8;X1]~[#7;X3](~[#8;X1])~[#8;X1]'))
        for match in match_NO3m:
            marked_as_visited += list(match)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetBondBetweenAtoms(match[1],match[0]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetAtomWithIdx(match[1]).SetFormalCharge(+1)
            self.mol.GetAtomWithIdx(match[0]).SetFormalCharge(-1)
            self.mol.GetAtomWithIdx(match[2]).SetFormalCharge(-1)
            charge -= 1

        match_SO4m = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#8;X1]~[#16;X4](~[#8;X1])(~[#8;X1])~[#8;X1]'))
        for match in match_SO4m:
            marked_as_visited += list(match)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[0]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetBondBetweenAtoms(match[1],match[4]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetAtomWithIdx(match[1]).SetFormalCharge(+1)
            self.mol.GetAtomWithIdx(match[0]).SetFormalCharge(-1)
            self.mol.GetAtomWithIdx(match[4]).SetFormalCharge(-1)
            charge -= 1

        match_SO3H = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[!#8]~[#16;X4](~[#8;X1])(~[#8;X1])~[#8;X2;H1]'))
        for match in match_SO3H:
            marked_as_visited += list(match)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[0]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetBondBetweenAtoms(match[1],match[4]).SetBondType(Chem.BondType.SINGLE)

        match_dmso = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]~[#16;X3](~[#8;X1])~[#6]'))
        for match in match_dmso:
            marked_as_visited += [match[1], match[2]]
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)

        match_dmso2 = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6,#7]~[#16;X4](~[#8;X1])(~[#8;X1])~[#6]'))
        for match in match_dmso2:
            marked_as_visited += [match[1], match[2], match[3]]
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.DOUBLE)

        match_cco = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[#6;X3;v3]-[#6;X3;v4]=[#8;X1;v2]'))
        first_idx = []
        for match in match_cco:
            if match[0] not in first_idx:
                first_idx.append(match[0])
                self.mol.GetBondBetweenAtoms(match[0],match[1]).SetBondType(Chem.BondType.DOUBLE)
                self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.SINGLE)
                self.mol.GetAtomWithIdx(match[2]).SetFormalCharge(-1)
                charge -= 1
                marked_as_visited += match

        match_cooh = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[!#8]~[#6;X3](~[#8;X1])~[#8;X2]~[#1;X1]'))
        for match in match_cooh:
            marked_as_visited += list(match)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetAtomWithIdx(match[1]).SetFormalCharge(0)
            self.mol.GetAtomWithIdx(match[2]).SetFormalCharge(0)
            self.mol.GetAtomWithIdx(match[3]).SetFormalCharge(0)

        match_coom = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[!#8]~[#6;X3](~[#8;X1])~[#8;X1]'))
        for match in match_coom:
            marked_as_visited += list(match)
            self.mol.GetBondBetweenAtoms(match[1],match[2]).SetBondType(Chem.BondType.DOUBLE)
            self.mol.GetBondBetweenAtoms(match[1],match[3]).SetBondType(Chem.BondType.SINGLE)
            self.mol.GetAtomWithIdx(match[1]).SetFormalCharge(0)
            self.mol.GetAtomWithIdx(match[2]).SetFormalCharge(0)
            self.mol.GetAtomWithIdx(match[3]).SetFormalCharge(-1)
            charge -= 1
        unvisited_charge = self.get_total_charge(marked_as_visited)
        self.charge = unvisited_charge + charge

        try:
            rdDetermineBonds.DetermineBonds(self.mol, charge=self.charge, useHueckel=False, useVdw=True)
            self.sanitized = True
        except ValueError as e:
            raise e

    def get_init_idxs(self, mol_idxs):
        return[self.idxs[mol_idx] for mol_idx in mol_idxs]

    def detect_group(self, idx, visited_sites, final_idxs, final_charge):
        if idx in visited_sites:
            return
        final_idxs.append(self.idxs[idx])
        a = self.mol.GetAtoms()[idx]
        final_charge[0] += a.GetFormalCharge()
        visited_sites.append(idx)
        for neighbor in a.GetNeighbors():
            self.detect_group(neighbor.GetIdx(), visited_sites, final_idxs, final_charge)

    def getect_grouped_charges(self):
        self.grouped_charges = []
        ### one_atom_ligand rdkit do not detect charges
        if not self.sanitized:
            self.grouped_charges.append(GroupedCharge(self.idxs, self.charge))
            return
        ### CO2 case
        matches_co2 = self.mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O-]"))
        for match in matches_co2:
            self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))

        matcher_adder = []
        ### NmNNN aromatic 5-member ring case
        matches_nm_nnn =  self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n-;R1;r5]~[n;R1;r5;X2]~[n;R1;r5;X2]~[n;R1;r5;X2]"))
        if len(matches_nm_nnn)>0:
            for match in matches_nm_nnn:
                self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
            matcher_adder += list(matches_nm_nnn)
        else:
            matches_nnm_nn =  self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n;R1;r5]~[n-;R1;r5;X2]~[n;R1;r5;X2]~[n;R1;r5;X2]"))
            if len(matches_nnm_nn)>0:
                for match in matches_nnm_nn:
                    self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
                matcher_adder += list(matches_nnm_nn)
            else:
                matches_nm_nn = self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n-;R1;r5]~[n;R1;r5;X2]~[n;R1;r5;X2]"))
                if len(matches_nm_nn)>0:
                    for match in matches_nm_nn:
                        self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
                    matcher_adder += list(matches_nm_nn)
                matches_nnm_n = self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n;R1;r5;X2]~[n-;R1;r5]~[n;R1;r5;X2]"))
                if len(matches_nnm_n)>0:
                    for match in matches_nnm_n:
                        self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
                    matcher_adder += list(matches_nnm_n)
                if len(matches_nnm_n)==0 and len(matches_nm_nn)==0:
                    matches_nm_cn =  self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n-;R1;r5][c;R1;r5][n;R1;r5;X2]"))
                    if len(matches_nm_cn)>0:
                        for match in matches_nm_cn:
                            self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
                        matcher_adder += list(matches_nm_cn)
                    matches_nn =  self.mol.GetSubstructMatches(Chem.MolFromSmarts("[n-;R1;r5][n;R1;r5;X2]"))
                    if len(matches_nn)>0:
                        for match in matches_nn:
                            self.grouped_charges.append(GroupedCharge(self.get_init_idxs(match), -1))
                        matcher_adder += list(matches_nn)

        visited_sites = [i for match in list(matches_co2) + matcher_adder for i in match]
        for atom in self.mol.GetAtoms():
            if atom.GetFormalCharge()==0:
                visited_sites.append(atom.GetIdx())
        for atom in self.mol.GetAtoms():
            final_idxs = []
            final_charge = [0.0]
            if atom.GetFormalCharge()!=0:
                self.detect_group(atom.GetIdx(), visited_sites, final_idxs, final_charge)
                self.grouped_charges.append(GroupedCharge(final_idxs, final_charge[0]))

    def detect_ions(self):
        for anion in Anions:
            if self.N == anion.N:
                match = self.mol.GetSubstructMatches(anion.mol)
                if len(match)>0:
                    if anion.total_charge != 0:
                        return True, GroupedCharge([self.idxs[i] for i in anion.charged_atoms], anion.total_charge)
                    else:
                        return True, None
        if self.N <= 3:
            raise Exception(f"Error. Uknown Anion {'_'.join([a.GetSymbol() for a in self.mol.GetAtoms()])}")
        else:
            return False, None


    def __init__(self, ligand_info_array, s, counter):
        logger.info(f"Working with ligand # {counter}")
        Zs = []
        coords = []
        symbols = []
        idxs = []
        for atom in ligand_info_array:
            idx = atom[0]
            Z = s[idx].specie.Z
            symbol = s[idx].specie.symbol
            frac_coords = s[idx].frac_coords + atom[1:]
            idxs.append(idx)
            Zs.append(Z)
            symbols.append(symbol)
            coords.append(frac_coords)
        coords = s.lattice.get_cartesian_coords(coords)
        self.idxs = idxs
        self.N = len(self.idxs)
        self.Zs = Zs
        self.symbols = symbols
        self.coords = coords
        self.gen_rdkit()
        logger.info(f"{counter} lig. rdkit Molecules Generated")
        try:
            adder, group = self.detect_ions()
            if adder:
                logger.info(f"{counter} lig. Seems to be known small molecule/anion")
                if group is not None:
                    self.grouped_charges = [group]
                else:
                    self.grouped_charges = []
            else:
                logger.info(f"{counter} lig. Seems to be uknown big molecule. Analyzing")
                self.sanitize()
                self.getect_grouped_charges()
            self.validity = True
        except Exception as e:
            logger.error(f"{counter} lig. {e}")
            self.validity = False
