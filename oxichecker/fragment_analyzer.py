# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
from pathlib import Path
from collections import defaultdict

FRAGMENT_CHARGES = {
   'N3': -1.0,
   'NO3': -1.0,
   'NO2': -1.0,
   'SCN': -1.0,
   'SeCN': -1.0,
   'CN': -1.0,
   'NCO': -1.0,
   'NH3': 0,
   'HCO2': -1.0,
   'CH3': -1.0,
   'C2H5': -1.0,
   'H2': 0.0,
   'O2': 0.0,
   'N2': 0.0
}

ANION_CHARGES = {
   'O': -2.0, 
   'F': -1.0, 
   'Cl': -1.0, 
   'Br': -1.0, 
   'I': -1.0, 
   'N': 3, 
   'S': -2,
   'Se': -2,
   'H': 1.0
}

METALS = {'Cd', 'Cu', 'Zn', 'Ni', 'Cr', 'V', 'Au', 'Na', 'Al', 'Mo', 'Fe', 'Dy', 'K', 'Co', 'Ag', 'Pr', 'Sc', 'Sr', 'Ho', 'Ca', 'Li', 'Zr', 'Gd', 'Pb', 'Ga', 'Rh', 'Lu', 'Tm', 'Nb', 'Th', 'Pu', 'Ir', 'Mn', 'Ba', 'Mg', 'Hg', 'Sn', 'U', 'La', 'Yb', 'Nd', 'Rb', 'Sm', 'Er', 'Eu', 'Tb', 'Tl', 'Ce', 'Ti', 'Pt', 'Ru', 'Bi', 'In', 'W', 'Cs', 'Pd', 'Hf', 'Y', 'B'}

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

def build_anion_graph(atoms, bonds):
   graph = defaultdict(list)
   
   for bond in bonds:
       atom1 = bond['atom1']
       atom2 = bond['atom2']
       
       if atoms[atom1]['element'] in METALS or atoms[atom2]['element'] in METALS:
           continue
       
       graph[atom1].append(atom2)
       graph[atom2].append(atom1)
   
   return graph

def find_anion_fragments(graph, atoms, max_iterations=10000):
   visited = set()
   fragments = []
   total_iterations = 0
   
   for atom_id in range(len(atoms)):
       if atom_id not in visited and atoms[atom_id]['element'] not in METALS:
           fragment = []
           stack = [atom_id]
           iterations = 0
           
           while stack and iterations < max_iterations:
               iterations += 1
               total_iterations += 1
               
               if total_iterations > max_iterations * 2:
                   print(f"Warning: Maximum iterations exceeded in find_anion_fragments")
                   break
               
               current = stack.pop()
               if current not in visited:
                   visited.add(current)
                   fragment.append(current)
                   # Only add unvisited neighbors to avoid cycles
                   for neighbor in graph[current]:
                       if neighbor not in visited:
                           stack.append(neighbor)
           
           if fragment:
               fragments.append(fragment)
   
   return fragments

def classify_fragment(fragment, atoms):
   composition = defaultdict(int)
   
   for atom_id in fragment:
       element = atoms[atom_id]['element']
       composition[element] += 1
   
   # Check for SeCN first (priority over CN)
   if composition.get('C', 0) == 1 and composition.get('N', 0) == 1 and composition.get('Se', 0) == 1:
       return 'SeCN'
   # Check for SCN (priority over CN)
   elif composition.get('C', 0) == 1 and composition.get('N', 0) == 1 and composition.get('S', 0) == 1:
       return 'SCN'
   # Check for NCO
   elif composition.get('N', 0) == 1 and composition.get('C', 0) == 1 and composition.get('O', 0) == 1:
       return 'NCO'
   # Check for CN (only if not part of SCN, SeCN or NCO)
   elif composition.get('C', 0) == 1 and composition.get('N', 0) == 1 and len(composition) == 2:
       return 'CN'
   # Check for NO3 (nitrate)
   elif composition.get('N', 0) == 1 and composition.get('O', 0) == 3 and len(composition) == 2:
       return 'NO3'
   # Check for NO2 (nitrite)
   elif composition.get('N', 0) == 1 and composition.get('O', 0) == 2 and len(composition) == 2:
       return 'NO2'
   # Check for N3 (azide) - must be checked before N2
   elif composition.get('N', 0) == 3 and len(composition) == 1:
       return 'N3'
   # Check for N2 (dinitrogen) - only if not N3
   elif composition.get('N', 0) == 2 and len(composition) == 1:
       return 'N2'
   # Check for O2 (dioxygen)
   elif composition.get('O', 0) == 2 and len(composition) == 1:
       return 'O2'
   # Check for H2 (dihydrogen)
   elif composition.get('H', 0) == 2 and len(composition) == 1:
       return 'H2'
   # Check for NH3 (ammonia)
   elif composition.get('N', 0) == 1 and composition.get('H', 0) == 3 and len(composition) == 2:
       return 'NH3'
   # Check for HCO2 (formate)
   elif composition.get('C', 0) == 1 and composition.get('H', 0) == 1 and composition.get('O', 0) == 2:
       return 'HCO2'
   # Check for CH3 (methyl)
   elif composition.get('C', 0) == 1 and composition.get('H', 0) == 3 and len(composition) == 2:
       return 'CH3'
   # Check for C2H5 (ethyl)
   elif composition.get('C', 0) == 2 and composition.get('H', 0) == 5 and len(composition) == 2:
       return 'C2H5'
   
   return None

def analyze_node_fragments(xyz_file):
   xyz_path = Path(xyz_file)
   
   if not xyz_path.exists():
       return 0.0, []
   
   atoms, bonds = parse_xyz(xyz_path)
   
   if not atoms:
       return 0.0, []
   
   graph = build_anion_graph(atoms, bonds)
   fragments = find_anion_fragments(graph, atoms)
   
   used_atoms = set()
   found_fragments = []
   total_charge = 0.0
   
   for fragment in fragments:
       ligand_type = classify_fragment(fragment, atoms)
       if ligand_type:
           fragment_charge = FRAGMENT_CHARGES[ligand_type]
           found_fragments.append({
               'name': ligand_type,
               'charge': fragment_charge,
               'atoms': [atoms[i]['element'] for i in fragment]
           })
           total_charge += fragment_charge
           used_atoms.update(fragment)
   
   for atom_id, atom in enumerate(atoms):
       if atom_id not in used_atoms and atom['element'] not in METALS:
           element = atom['element']
           if element in ANION_CHARGES:
               atom_charge = ANION_CHARGES[element]
               found_fragments.append({
                   'name': f'isolated_{element}',
                   'charge': atom_charge,
                   'atoms': [element]
               })
               total_charge += atom_charge
   
   return total_charge, found_fragments