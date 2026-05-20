# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
from pathlib import Path
from collections import defaultdict

PILLAR_METALS = {'Nb', 'Ge', 'Ti', 'Si', 'Zr', 'Fe', 'Al', 'V', 'Sn', 'Ga', 'In', 'B', 'Ir'}
PILLAR_LIGANDS = {'O', 'F'}
PILLAR_CHARGE = -2.0

PERCHLORATE_CENTERS = {'Cl', 'Br'}
PERCHLORATE_CHARGE = -1.0

def parse_xyz(xyz_file):
   atoms = []
   bonds = []
   
   with open(xyz_file, 'r') as f:
       lines = f.readlines()
   
   n_atoms = int(lines[0].strip())
   
   for i in range(2, 2 + n_atoms):
       parts = lines[i].split()
       element = parts[0]
       atoms.append({'id': len(atoms), 'element': element})
   
   for i in range(2 + n_atoms, len(lines)):
       parts = lines[i].split()
       if len(parts) >= 2:
           atom1 = int(parts[0])
           atom2 = int(parts[1])
           bonds.append({'atom1': atom1, 'atom2': atom2})
   
   return atoms, bonds

def build_connectivity_graph(atoms, bonds):
   graph = defaultdict(list)
   
   for bond in bonds:
       atom1 = bond['atom1']
       atom2 = bond['atom2']
       graph[atom1].append(atom2)
       graph[atom2].append(atom1)
   
   return graph

def check_pillar_complex(atoms, graph):
   for atom_id, atom in enumerate(atoms):
       if atom['element'] in PILLAR_METALS:
           neighbors = list(set(graph[atom_id]))
           
           neighbor_composition = defaultdict(int)
           for neighbor_id in neighbors:
               neighbor_element = atoms[neighbor_id]['element']
               if neighbor_element in PILLAR_LIGANDS:
                   neighbor_composition[neighbor_element] += 1
           
           num_O = neighbor_composition.get('O', 0)
           num_F = neighbor_composition.get('F', 0)
           total_ligands = num_O + num_F
           
           if (total_ligands == 6 and num_F >= 1) or (total_ligands == 5 and num_F >= 1):
               return ('replace', PILLAR_CHARGE)
   
   return None

def check_perchlorate(atoms, graph):
   for atom_id, atom in enumerate(atoms):
       if atom['element'] in PERCHLORATE_CENTERS:
           neighbors = list(set(graph[atom_id]))
           
           o_count = 0
           for neighbor_id in neighbors:
               if atoms[neighbor_id]['element'] == 'O':
                   o_count += 1
           
           if o_count == 4:
               return ('replace', PERCHLORATE_CHARGE)
   
   return None

def count_PO3C_groups(atoms, graph):
   po3c_count = 0
   
   for atom_id, atom in enumerate(atoms):
       if atom['element'] == 'P':
           neighbors = list(set(graph[atom_id]))
           
           o_neighbors = []
           has_carbon = False
           
           for neighbor_id in neighbors:
               neighbor_element = atoms[neighbor_id]['element']
               if neighbor_element == 'O':
                   o_neighbors.append(neighbor_id)
               elif neighbor_element == 'C':
                   has_carbon = True
           
           if len(o_neighbors) == 3 and has_carbon:
               has_oh_group = False
               for o_id in o_neighbors:
                   o_neighbors_of_o = graph[o_id]
                   for neighbor_of_o in o_neighbors_of_o:
                       if atoms[neighbor_of_o]['element'] == 'H':
                           has_oh_group = True
                           break
                   if has_oh_group:
                       break
               
               if not has_oh_group:
                   po3c_count += 1
   
   return po3c_count

def count_SO3_groups(atoms, graph):
   so3_count = 0
   
   for atom_id, atom in enumerate(atoms):
       if atom['element'] == 'S':
           neighbors = list(set(graph[atom_id]))
           
           has_s_neighbor = any(atoms[n_id]['element'] == 'S' for n_id in neighbors)
           if has_s_neighbor:
               continue
           
           o_neighbors = []
           for n_id in neighbors:
               if atoms[n_id]['element'] == 'O':
                   o_neighbors.append(n_id)
           
           if len(o_neighbors) == 3:
               has_oh_group = False
               for o_id in o_neighbors:
                   o_neighbors_of_o = graph[o_id]
                   for neighbor_of_o in o_neighbors_of_o:
                       if atoms[neighbor_of_o]['element'] == 'H':
                           has_oh_group = True
                           break
                   if has_oh_group:
                       break
               
               if not has_oh_group:
                   so3_count += 1
   
   return so3_count

def has_PO3_groups(xyz_file):
   xyz_path = Path(xyz_file)
   
   if not xyz_path.exists():
       return False
   
   atoms, bonds = parse_xyz(xyz_path)
   
   if not atoms:
       return False
   
   graph = build_connectivity_graph(atoms, bonds)
   po3c_count = count_PO3C_groups(atoms, graph)
   
   return po3c_count > 0

def analyze_linker_charge(xyz_file):
   xyz_path = Path(xyz_file)
   
   if not xyz_path.exists():
       return None
   
   atoms, bonds = parse_xyz(xyz_path)
   
   if not atoms:
       return None
   
   graph = build_connectivity_graph(atoms, bonds)
   
   result = check_pillar_complex(atoms, graph)
   if result is not None:
       return result
   
   result = check_perchlorate(atoms, graph)
   if result is not None:
       return result
   
   po3c_count = count_PO3C_groups(atoms, graph)
   so3_count = count_SO3_groups(atoms, graph)
   
   total_adjustment = 0
   
#    if po3c_count > 0:
#        total_adjustment -= po3c_count * 1.0
   
   if so3_count > 0:
       total_adjustment += so3_count * 2
   
   if total_adjustment != 0:
       return ('adjust', total_adjustment)
   
   return None