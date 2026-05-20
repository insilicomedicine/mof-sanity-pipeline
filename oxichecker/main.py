#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import subprocess
from pathlib import Path
import sys
import argparse
import csv
from itertools import product
from rdkit import Chem
from rdkit.Chem import rdmolops
import tempfile
import shutil
from multiprocessing import Process, Queue, cpu_count
from tqdm import tqdm
from functools import partial
import warnings
from rdkit import RDLogger
import gc
import time
from fragment_analyzer import analyze_node_fragments
from linker_analyzer import analyze_linker_charge, has_PO3_groups

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=SyntaxWarning)
RDLogger.DisableLog('rdApp.*')


CIF_TO_BB_PATH = Path("/opt/libcif/bin/cif_to_building_blocks_v3")

def is_node_block(block_name):
    name_lower = block_name.lower()
    
    if "linker" in name_lower or "disconnected" in name_lower:
        return False
    
    if "node" in name_lower:
        return True
    
    return False


def is_linker_block(block_name):
    name_lower = block_name.lower()
    return "linker" in name_lower or "disconnected" in name_lower

def has_metals_in_cif(cif_path):
    """Check if CIF file contains any metals from NORMAL_CHARGES"""
    try:
        with open(cif_path, 'r') as f:
            for line in f:
                line = line.strip()
                # Look for _atom_site_type_symbol or similar atom declarations
                if line.startswith('_atom_site_type_symbol') or line.startswith('_chemical_formula'):
                    continue
                # Check data lines that might contain element symbols
                if line and not line.startswith('_') and not line.startswith('#'):
                    parts = line.split()
                    if parts:
                        # First part might be element symbol
                        element = parts[0].strip()
                        # Remove numbers and special characters
                        clean_element = ''.join(c for c in element if c.isalpha())
                        if clean_element in NORMAL_CHARGES:
                            return True
        return False
    except:
        return True  # If can't read file, assume it has metals to be safe

NORMAL_CHARGES = {
   'Cu': [1.0, 2.0], 'Zn': [2.0], 'Ni': [2.0, 3.0], 'Cr': [3.0, 6.0], 'V': [3.0, 5.0], 'Au': [1.0, 3.0],
   'Na': [1.0], 'Cd': [2.0, 3.0], 'Al': [3.0], 'Mo': [4.0, 6.0], 'Fe': [2.0, 3.0], 'Dy': [3.0],
   'K': [1.0], 'Co': [2.0, 3.0], 'Ag': [1.0], 'Pr': [3.0], 'Sc': [3.0], 'Sr': [2.0], 'Ho': [3.0],
   'Ca': [2.0], 'Li': [1.0], 'Zr': [2.0, 3.0, 4.0], 'Gd': [3.0], 'Pb': [2.0, 4.0], 'Ga': [3.0], 'Rh': [3.0], 'Ge': [2.0, 4.0],
   'Lu': [3.0], 'Tm': [3.0], 'Nb': [5.0], 'Th': [4.0], 'Pu': [4.0], 'Ir': [3.0, 4.0], 'Sb': [3.0, 5.0],
   'Mn': [2.0, 3.0], 'Ba': [2.0], 'Mg': [2.0], 'Hg': [2.0], 'Sn': [0.0, 2.0, 4.0], 'U': [6.0],'Re': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
   'La': [3.0], 'Yb': [3.0], 'Nd': [3.0], 'Rb': [1.0], 'Sm': [3.0], 'Er': [3.0], 'Eu': [3.0],'Np': [3.0, 4.0, 5.0, 6.0, 7.0],
   'Tb': [3.0], 'Tl': [3.0], 'Ce': [1.0], 'Ti': [2.0, 4.0], 'Pt': [2.0, 4.0], 'Ru': [3.0, 4.0],
   'Bi': [3.0], 'In': [3.0], 'W': [6.0], 'Cs': [1.0], 'Pd': [2.0], 'Hf': [4.0], 'Y': [3.0], "B": [0.0, 3.0]
}


def run_cif_decomposition(cif_path, work_dir, cif_to_bb_path, decomp_params, verbose=False, timeout=60):
    cif_filename = Path(cif_path).name
    cmd = [str(cif_to_bb_path), cif_filename] + decomp_params
    
    if verbose:
        print(f"[VERBOSE] Running: {' '.join(cmd)}")
        print(f"[VERBOSE] Working directory: {work_dir}")
        print(f"[VERBOSE] CIF file: {cif_path}")
    
    try:
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=timeout)
        
        if verbose:
            print(f"[VERBOSE] Exit code: {result.returncode}")
            if result.stdout:
                print(f"[VERBOSE] STDOUT: {result.stdout}")
            if result.stderr:
                print(f"[VERBOSE] STDERR: {result.stderr}")
        
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        if verbose:
            print(f"[VERBOSE] TIMEOUT: Process killed after 60 seconds for {cif_path}")
        return False
    except Exception as e:
        if verbose:
            print(f"[VERBOSE] ERROR: {e} for {cif_path}")
        return False

def read_building_blocks(work_dir):
   num_blocks_file = work_dir / "num_building_blocks.txt"
   if not num_blocks_file.exists():
       return None
   
   blocks = {}
   with open(num_blocks_file, 'r') as f:
       for line in f:
           if ':' in line:
               name, count = line.strip().split(':')
               blocks[name.strip()] = int(count.strip())
   
   return blocks


def read_infinite_node_blocks(work_dir, node_name):
   infinite_blocks_file = work_dir / f"{node_name}_num_building_blocks.txt"
   if not infinite_blocks_file.exists():
       return None
   
   blocks = {}
   with open(infinite_blocks_file, 'r') as f:
       for line in f:
           if ':' in line:
               name, count = line.strip().split(':')
               blocks[name.strip()] = int(count.strip())
   
   return blocks


def process_infinite_nodes(temp_path, blocks):
   expanded_blocks = {}
   infinite_nodes_to_skip = set()
   
   for block_name, count in blocks.items():
       if "_infinite" in block_name and is_node_block(block_name):
           node_base = block_name.replace(".xyz", "")
           
           infinite_node_blocks = read_infinite_node_blocks(temp_path, node_base)
           
           if infinite_node_blocks:
               infinite_nodes_to_skip.add(block_name)
               
               for sub_block_name, sub_count in infinite_node_blocks.items():
                   full_name = f"{node_base}_{sub_block_name}"
                   total_count = sub_count * count
                   
                   if full_name in expanded_blocks:
                       expanded_blocks[full_name] += total_count
                   else:
                       expanded_blocks[full_name] = total_count
           else:
               expanded_blocks[block_name] = count
       else:
           expanded_blocks[block_name] = count
   
   return expanded_blocks, infinite_nodes_to_skip


def get_node_elements(xyz_file):
   elements = []
   with open(xyz_file, 'r') as f:
       lines = f.readlines()
   
   n_atoms = int(lines[0].strip())
   for i in range(2, min(2 + n_atoms, len(lines))):
       parts = lines[i].split()
       if parts:
           elements.append(parts[0])
   
   return elements


def path1_xyz_to_inchi_to_rdkit(xyz_file, timeout=30):
    try:
        result = subprocess.run(['obabel', '-ixyz', str(xyz_file), '-oinchi', '-xu', '-as'], 
                              capture_output=True, text=True, timeout=timeout)
        if result.returncode == 0:
            inchi = result.stdout.strip()
            mol = Chem.MolFromInchi(inchi)
            return mol, inchi
    except subprocess.TimeoutExpired:
        pass
    except Exception:
        pass
    return None, ""


def path2_xyz_to_mol2_to_inchi_to_rdkit(xyz_file, timeout=30):
    mol2_file = xyz_file.with_suffix('.mol2')
    
    try:
        cmd_mol2 = ["obabel", "-ixyz", str(xyz_file), "-omol2", "-O", str(mol2_file)]
        result = subprocess.run(cmd_mol2, capture_output=True, timeout=timeout)
        
        if not mol2_file.exists():
            return None, ""
        
        result = subprocess.run(['obabel', '-imol2', str(mol2_file), '-oinchi'], 
                              capture_output=True, text=True, timeout=timeout)
        
        mol2_file.unlink()
        
        if result.returncode == 0:
            inchi = result.stdout.strip()
            mol = Chem.MolFromInchi(inchi)
            return mol, inchi
    except subprocess.TimeoutExpired:
        if mol2_file.exists():
            mol2_file.unlink()
    except Exception:
        if mol2_file.exists():
            mol2_file.unlink()
    
    return None, ""


def path3_xyz_to_mol_to_rdkit_to_smiles(xyz_file, timeout=30):
    mol_file = xyz_file.with_suffix('.mol')
    
    try:
        cmd_mol = ["obabel", "-ixyz", str(xyz_file), "-omol", "-O", str(mol_file)]
        result = subprocess.run(cmd_mol, capture_output=True, timeout=timeout)
        
        if not mol_file.exists():
            return None, ""
        
        mol = Chem.MolFromMolFile(str(mol_file))
        mol_file.unlink()
        
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            return mol, smiles
    except subprocess.TimeoutExpired:
        if mol_file.exists():
            mol_file.unlink()
    except Exception:
        if mol_file.exists():
            mol_file.unlink()
    
    return None, ""


def calculate_charge_from_rdkit_mol(mol):
   if mol is None:
       return 0
   
   formal_charge = rdmolops.GetFormalCharge(mol)
   radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
   total_charge = formal_charge - radicals
   
   return total_charge

def validate_mof_worker(args_tuple):
    cif_path, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout = args_tuple
    
    # Clean CIF name by removing P1_ prefix if present
    cif_name = cif_path.name
    if cif_name.startswith('P1_'):
        cif_name = cif_name[3:]  # Remove 'P1_' prefix


    if decompose_mode:
        structure_name = cif_path.stem
        temp_path = output_base_dir / structure_name
        temp_path.mkdir(exist_ok=True, mode=0o777)
        temp_cif = temp_path / cif_path.name
        shutil.copy2(cif_path, temp_cif)
        cleanup_temp = False
    else:
        temp_dir = tempfile.TemporaryDirectory()
        temp_path = Path(temp_dir.name)
        temp_cif = temp_path / cif_path.name
        shutil.copy2(cif_path, temp_cif)
        cleanup_temp = True
        
    try:
        # Check if CIF contains metals before attempting decomposition
        if not has_metals_in_cif(temp_cif):
            return {
                'cif': cif_name,
                'OxiChecker Validity': 'True; no metal',
                'PATH1_InChI': '',
                'PATH2_InChI': '',
                'PATH3_SMILES': '',
                'Valid_Path': ''
            }
        
        if not run_cif_decomposition(temp_cif, temp_path, cif_to_bb_path, decomp_params, verbose, decomp_timeout):
            return {
                'cif': cif_name,
                'OxiChecker Validity': 'Failed to decompose',
                'PATH1_InChI': '',
                'PATH2_InChI': '',
                'PATH3_SMILES': '',
                'Valid_Path': ''
            }
        
        blocks = read_building_blocks(temp_path)
        if not blocks:
            xyz_files = list(temp_path.glob("*.xyz"))
            if not xyz_files:
                return {
                    'cif': cif_name,
                    'OxiChecker Validity': 'Empty decomposition',
                    'PATH1_InChI': '',
                    'PATH2_InChI': '',
                    'PATH3_SMILES': '',
                    'Valid_Path': ''
                }
            else:
                return {
                    'cif': cif_name,
                    'OxiChecker Validity': 'No building blocks file',
                    'PATH1_InChI': '',
                    'PATH2_InChI': '',
                    'PATH3_SMILES': '',
                    'Valid_Path': ''
                }
        
        if all(count == 0 for count in blocks.values()):
            return {
                'cif': cif_name,
                'OxiChecker Validity': 'Empty decomposition',
                'PATH1_InChI': '',
                'PATH2_InChI': '',
                'PATH3_SMILES': '',
                'Valid_Path': ''
            }
        
        expanded_blocks, infinite_nodes_to_skip = process_infinite_nodes(temp_path, blocks)
        
        infinite_nodes = []
        infinite_linkers = []
        
        for block_name in expanded_blocks.keys():
            if "_infinite" in block_name and block_name not in infinite_nodes_to_skip:
                if is_node_block(block_name):
                    infinite_nodes.append(block_name)
                elif is_linker_block(block_name):
                    infinite_linkers.append(block_name)
        
        linker_data_path1 = []
        linker_data_path2 = []
        linker_data_path3 = []
        node_data = []
        
        linker_inchi_path1 = []
        linker_inchi_path2 = []
        linker_smiles_path3 = []
        
        for block_name, count in expanded_blocks.items():
            block_file = temp_path / block_name
            
            if is_linker_block(block_name) and block_file.exists():
                if has_PO3_groups(block_file):
                    mol3, smiles3 = path3_xyz_to_mol_to_rdkit_to_smiles(block_file, obabel_timeout)
                    
                    if not smiles3 or not mol3:
                        mol1, inchi1 = path1_xyz_to_inchi_to_rdkit(block_file, obabel_timeout)
                        mol2, inchi2 = path2_xyz_to_mol2_to_inchi_to_rdkit(block_file, obabel_timeout)
                        
                        linker_inchi_path1.append(f"{block_name}: {inchi1}" if inchi1 else f"{block_name}: None")
                        linker_inchi_path2.append(f"{block_name}: {inchi2}" if inchi2 else f"{block_name}: None")
                        linker_smiles_path3.append(f"{block_name}: None")
                        
                        special_result = analyze_linker_charge(block_file)
                        
                        if special_result and special_result[0] == 'replace':
                            charge_value = special_result[1]
                            if mol1:
                                linker_data_path1.append({'charge': charge_value, 'count': count, 'name': block_name})
                            if mol2:
                                linker_data_path2.append({'charge': charge_value, 'count': count, 'name': block_name})
                        elif special_result and special_result[0] == 'adjust':
                            adjustment = special_result[1]
                            if mol1:
                                charge1 = calculate_charge_from_rdkit_mol(mol1) + adjustment
                                linker_data_path1.append({'charge': charge1, 'count': count, 'name': block_name})
                            if mol2:
                                charge2 = calculate_charge_from_rdkit_mol(mol2) + adjustment
                                linker_data_path2.append({'charge': charge2, 'count': count, 'name': block_name})
                        else:
                            if mol1:
                                charge1 = calculate_charge_from_rdkit_mol(mol1)
                                linker_data_path1.append({'charge': charge1, 'count': count, 'name': block_name})
                            if mol2:
                                charge2 = calculate_charge_from_rdkit_mol(mol2)
                                linker_data_path2.append({'charge': charge2, 'count': count, 'name': block_name})
                    else:
                        # PATH3 succeeded for PO3
                        linker_smiles_path3.append(f"{block_name}: {smiles3}")
                        linker_inchi_path1.append(f"{block_name}: PO3_group_PATH3_only")
                        linker_inchi_path2.append(f"{block_name}: PO3_group_PATH3_only")
                        
                        base_charge3 = calculate_charge_from_rdkit_mol(mol3)
                        
                        special_result = analyze_linker_charge(block_file)
                        final_charge3 = base_charge3
                        
                        if special_result and special_result[0] == 'adjust':
                            adjustment = special_result[1]
                            final_charge3 = base_charge3 + adjustment
                        
                        linker_data_path3.append({'charge': final_charge3, 'count': count, 'name': block_name})
                    continue
                
                special_result = analyze_linker_charge(block_file)
                if special_result is not None:
                    action, value = special_result
                    
                    if action == 'replace':
                        linker_data_path1.append({'charge': value, 'count': count, 'name': block_name})
                        linker_data_path2.append({'charge': value, 'count': count, 'name': block_name})
                        linker_data_path3.append({'charge': value, 'count': count, 'name': block_name})
                        
                        mol1, inchi1 = path1_xyz_to_inchi_to_rdkit(block_file, obabel_timeout)
                        mol2, inchi2 = path2_xyz_to_mol2_to_inchi_to_rdkit(block_file, obabel_timeout)
                        mol3, smiles3 = path3_xyz_to_mol_to_rdkit_to_smiles(block_file, obabel_timeout)
                        
                        linker_inchi_path1.append(f"{block_name}: {inchi1}" if inchi1 else f"{block_name}: None")
                        linker_inchi_path2.append(f"{block_name}: {inchi2}" if inchi2 else f"{block_name}: None")
                        linker_smiles_path3.append(f"{block_name}: {smiles3}" if smiles3 else f"{block_name}: None")
                        continue
                    
                    elif action == 'adjust':
                        adjustment = value
                        
                        mol1, inchi1 = path1_xyz_to_inchi_to_rdkit(block_file, obabel_timeout)
                        if mol1:
                            charge1 = calculate_charge_from_rdkit_mol(mol1) + adjustment
                            linker_data_path1.append({'charge': charge1, 'count': count, 'name': block_name})
                        
                        mol2, inchi2 = path2_xyz_to_mol2_to_inchi_to_rdkit(block_file, obabel_timeout)
                        if mol2:
                            charge2 = calculate_charge_from_rdkit_mol(mol2) + adjustment
                            linker_data_path2.append({'charge': charge2, 'count': count, 'name': block_name})
                        
                        mol3, smiles3 = path3_xyz_to_mol_to_rdkit_to_smiles(block_file, obabel_timeout)
                        if mol3:
                            charge3 = calculate_charge_from_rdkit_mol(mol3) + adjustment
                            linker_data_path3.append({'charge': charge3, 'count': count, 'name': block_name})
                        
                        linker_inchi_path1.append(f"{block_name}: {inchi1}" if inchi1 else f"{block_name}: None")
                        linker_inchi_path2.append(f"{block_name}: {inchi2}" if inchi2 else f"{block_name}: None")
                        linker_smiles_path3.append(f"{block_name}: {smiles3}" if smiles3 else f"{block_name}: None")
                        
                        continue
                
                mol1, inchi1 = path1_xyz_to_inchi_to_rdkit(block_file, obabel_timeout)
                mol2, inchi2 = path2_xyz_to_mol2_to_inchi_to_rdkit(block_file, obabel_timeout)
                mol3, smiles3 = path3_xyz_to_mol_to_rdkit_to_smiles(block_file, obabel_timeout)
                
                linker_inchi_path1.append(f"{block_name}: {inchi1}" if inchi1 else f"{block_name}: None")
                linker_inchi_path2.append(f"{block_name}: {inchi2}" if inchi2 else f"{block_name}: None")
                linker_smiles_path3.append(f"{block_name}: {smiles3}" if smiles3 else f"{block_name}: None")
                
                if mol1:
                    charge1 = calculate_charge_from_rdkit_mol(mol1)
                    linker_data_path1.append({'charge': charge1, 'count': count, 'name': block_name})
                
                if mol2:
                    charge2 = calculate_charge_from_rdkit_mol(mol2)
                    linker_data_path2.append({'charge': charge2, 'count': count, 'name': block_name})
                
                if mol3:
                    charge3 = calculate_charge_from_rdkit_mol(mol3)
                    linker_data_path3.append({'charge': charge3, 'count': count, 'name': block_name})
        
        for block_name, count in expanded_blocks.items():
            block_file = temp_path / block_name
            
            if is_node_block(block_name) and block_file.exists():
                if block_name in infinite_nodes_to_skip:
                    continue
                    
                elements = get_node_elements(block_file)
                composition = {}
                for elem in elements:
                    composition[elem] = composition.get(elem, 0) + 1
                
                metals = {}
                for elem, elem_count in composition.items():
                    if elem in NORMAL_CHARGES:
                        metals[elem] = elem_count
                
                anion_charge, fragments = analyze_node_fragments(block_file)
                
                node_data.append({
                    'count': count,
                    'metals': metals,
                    'anion_charge': anion_charge,
                    'name': block_name
                })
        
        def check_path_validity(linker_data, path_name):
            total_linker_charge = sum(l['charge'] * l['count'] for l in linker_data)
            total_anion_charge = sum(n['anion_charge'] * n['count'] for n in node_data)
            required_metal_charge = -total_linker_charge - total_anion_charge
            
            all_metals = []
            for node in node_data:
                for metal, count in node['metals'].items():
                    all_metals.extend([metal] * count * node['count'])
            
            # Check if structure has no metals (e.g., COF)
            if not all_metals and abs(required_metal_charge) < 0.01:
                return 'no_metal', []
            
            if not node_data:
                return False, []
            
            all_node_options = []
            for node in node_data:
                node_options = []
                for metal, metal_count in node['metals'].items():
                    charges = NORMAL_CHARGES.get(metal, [])
                    node_options.append([(metal, charge, metal_count) for charge in charges])
                
                if node_options:
                    node_combos = []
                    for combo in product(*node_options):
                        total_charge_per_node = sum(charge * count for _, charge, count in combo)
                        total_charge_all_nodes = total_charge_per_node * node['count']
                        node_combos.append((total_charge_all_nodes, combo))
                    all_node_options.append(node_combos)
            
            if not all_node_options:
                return False, all_metals
            
            if len(all_node_options) == 1:
                for charge, combo in all_node_options[0]:
                    if abs(charge - required_metal_charge) < 0.01:
                        valid_charges = []
                        for metal, charge_val, _ in combo:
                            valid_charges.append(f"{metal}{'+' if charge_val > 0 else ''}{int(charge_val)}")
                        return True, valid_charges
            else:
                for combo in product(*all_node_options):
                    total_charge = sum(option[0] for option in combo)
                    if abs(total_charge - required_metal_charge) < 0.01:
                        valid_charges = []
                        for node_option in combo:
                            for metal, charge_val, _ in node_option[1]:
                                valid_charges.append(f"{metal}{'+' if charge_val > 0 else ''}{int(charge_val)}")
                        return True, valid_charges
            
            return False, all_metals

        def check_mixed_path_validity():
            linker_path_combos = []
            linker_has_po3 = {}
            for block_name, count in expanded_blocks.items():
                if is_linker_block(block_name):
                    linker_path_combos.append(block_name)
                    block_file = temp_path / block_name
                    if block_file.exists():
                        linker_has_po3[block_name] = has_PO3_groups(block_file)
                    else:
                        linker_has_po3[block_name] = False
            
            if not linker_path_combos:
                return False, [], ""
            
            for combo in product([1, 2, 3], repeat=len(linker_path_combos)):
                mixed_linker_data = []
                combo_name = []
                
                for i, block_name in enumerate(linker_path_combos):
                    if linker_has_po3.get(block_name, False):
                        path3_exists = False
                        for linker in linker_data_path3:
                            if linker['name'] == block_name:
                                path3_exists = True
                                break
                        
                        if path3_exists:
                            path_choice = 3
                        else:
                            path_choice = combo[i]
                    else:
                        path_choice = combo[i]
                    combo_name.append(f"{block_name}:PATH{path_choice}")
                    
                    if path_choice == 1:
                        for linker in linker_data_path1:
                            if linker['name'] == block_name:
                                mixed_linker_data.append(linker)
                                break
                    elif path_choice == 2:
                        for linker in linker_data_path2:
                            if linker['name'] == block_name:
                                mixed_linker_data.append(linker)
                                break
                    else:
                        for linker in linker_data_path3:
                            if linker['name'] == block_name:
                                mixed_linker_data.append(linker)
                                break
                
                if len(mixed_linker_data) == len(linker_path_combos):
                    total_linker_charge = sum(l['charge'] * l['count'] for l in mixed_linker_data)
                    total_anion_charge = sum(n['anion_charge'] * n['count'] for n in node_data)
                    required_metal_charge = -total_linker_charge - total_anion_charge
                    
                    all_metals = []
                    for node in node_data:
                        for metal, count in node['metals'].items():
                            all_metals.extend([metal] * count * node['count'])
                    
                    if node_data:
                        valid_combinations = []
                        for node in node_data:
                            node_options = []
                            for metal, metal_count in node['metals'].items():
                                charges = NORMAL_CHARGES.get(metal, [])
                                node_options.append([(metal, charge, metal_count) for charge in charges])
                            
                            if node_options:
                                for node_combo in product(*node_options):
                                    total_charge_per_node = sum(charge * count for _, charge, count in node_combo)
                                    total_charge_all_nodes = total_charge_per_node * node['count']
                                    valid_combinations.append((total_charge_all_nodes, node_combo, node['count']))
                        
                        if valid_combinations:
                            actual_node_count = len(node_data)
                            combo_count = len(valid_combinations)

                            if combo_count == actual_node_count:
                                if actual_node_count == 1:
                                    for combo_data in valid_combinations:
                                        total_charge_all_nodes, node_combo, node_count = combo_data
                                        if abs(total_charge_all_nodes - required_metal_charge) < 0.01:
                                            valid_charges = []
                                            for metal, charge, _ in node_combo:
                                                valid_charges.append(f"{metal}{'+' if charge > 0 else ''}{int(charge)}")
                                            return True, valid_charges, "+".join(combo_name)
                                else:
                                    total_available_charge = sum(combo_data[0] for combo_data in valid_combinations)
                                    
                                    if abs(total_available_charge - required_metal_charge) < 0.01:
                                        valid_charges = []
                                        for combo_data in valid_combinations:
                                            total_charge, node_combo, node_count = combo_data
                                            for metal, charge, _ in node_combo:
                                                valid_charges.append(f"{metal}{'+' if charge > 0 else ''}{int(charge)}")
                                        return True, valid_charges, "+".join(combo_name)
                            else:
                                for combo_data in valid_combinations:
                                    total_charge_all_nodes, node_combo, node_count = combo_data
                                    if abs(total_charge_all_nodes - required_metal_charge) < 0.01:
                                        valid_charges = []
                                        for metal, charge, _ in node_combo:
                                            valid_charges.append(f"{metal}{'+' if charge > 0 else ''}{int(charge)}")
                                        return True, valid_charges, "+".join(combo_name)
            
            return False, [], ""
        
        if not linker_data_path1 and not linker_data_path2 and not linker_data_path3:
            path1_result, path1_charges = check_path_validity([], "NO_LINKERS")
            path1_valid = path1_result == True or path1_result == 'no_metal'
            path2_valid, path2_charges = False, []
            path3_valid, path3_charges = False, []
            has_no_metal = path1_result == 'no_metal'
            mixed_valid, mixed_charges, mixed_path = False, [], ""
        else:
            path1_result, path1_charges = check_path_validity(linker_data_path1, "PATH1") if linker_data_path1 else (False, [])
            path2_result, path2_charges = check_path_validity(linker_data_path2, "PATH2") if linker_data_path2 else (False, [])
            path3_result, path3_charges = check_path_validity(linker_data_path3, "PATH3") if linker_data_path3 else (False, [])
            
            # Handle no_metal case
            path1_valid = path1_result == True or path1_result == 'no_metal'
            path2_valid = path2_result == True or path2_result == 'no_metal'
            path3_valid = path3_result == True or path3_result == 'no_metal'
            has_no_metal = path1_result == 'no_metal' or path2_result == 'no_metal' or path3_result == 'no_metal'
            
            mixed_valid, mixed_charges, mixed_path = check_mixed_path_validity()
        
        final_valid = path1_valid or path2_valid or path3_valid or mixed_valid
        
        valid_path_parts = []
        if linker_inchi_path1: 
            linker_names = [name.split(':')[0] for name in linker_inchi_path1]
            for linker_name in linker_names:
                block_file = temp_path / linker_name
                has_po3 = block_file.exists() and has_PO3_groups(block_file)
                
                if has_po3:
                    if path3_valid:
                        valid_paths = ["PATH3_PO3"]
                        alternative_paths = []
                    else:
                        valid_paths = []
                        alternative_paths = ["PATH3_PO3"]
                else:
                    valid_paths = []
                    alternative_paths = []
                    
                    if path1_valid:
                        valid_paths.append("PATH1")
                    else:
                        alternative_paths.append("PATH1")
                        
                    if path2_valid:
                        if not valid_paths:  
                            valid_paths.append("PATH2")
                        else:
                            alternative_paths.append("PATH2")
                    else:
                        alternative_paths.append("PATH2")
                        
                    if path3_valid:
                        if not valid_paths: 
                            valid_paths.append("PATH3")
                        else:
                            alternative_paths.append("PATH3")
                    else:
                        alternative_paths.append("PATH3")
                
                if mixed_valid and mixed_path:
                    valid_path_parts.append(f"{linker_name}:\{mixed_path}\\")
                else:
                    main_path = valid_paths[0] if valid_paths else "NONE"
                    alt_paths = ",".join(alternative_paths) if alternative_paths else ""
                    if alt_paths:
                        valid_path_parts.append(f"{linker_name}:\{main_path}({alt_paths})\\")
                    else:
                        valid_path_parts.append(f"{linker_name}:\{main_path}\\")
        
        valid_path_result = "; ".join(valid_path_parts)
        
        infinity_suffix = ""
        if infinite_nodes:
            infinity_suffix = "; Infinity Node"
        elif infinite_linkers:
            infinity_suffix = "; Infinity Linker"
        
        # Check if structure has no metals
        if 'has_no_metal' in locals() and has_no_metal:
            validity_result = 'True; no metal' + infinity_suffix
        else:
            validity_result = ('True' if final_valid else 'Invalid charges') + infinity_suffix

        return {
            'cif': cif_name,
            'OxiChecker Validity': validity_result,
            'PATH1_InChI': "; ".join(linker_inchi_path1) if linker_inchi_path1 else "",
            'PATH2_InChI': "; ".join(linker_inchi_path2) if linker_inchi_path2 else "",
            'PATH3_SMILES': "; ".join(linker_smiles_path3) if linker_smiles_path3 else "",
            'Valid_Path': valid_path_result
        }
    
    finally:
        if cleanup_temp:
            temp_dir.cleanup()


def worker_process(cif_file, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout, result_queue):
    # Clean CIF name by removing P1_ prefix if present
    cif_name = Path(cif_file).name 
    if cif_name.startswith('P1_'):
        cif_name = cif_name[3:]  # Remove 'P1_' prefix
    try:
        args_tuple = (cif_file, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout)
        result = validate_mof_worker(args_tuple)
        result_queue.put(result)
    except Exception as e:
        result = {
            'cif': cif_name,
            'OxiChecker Validity': f'Processing error: {str(e)}',
            'PATH1_InChI': '',
            'PATH2_InChI': '',
            'PATH3_SMILES': '',
            'Valid_Path': ''
        }
        result_queue.put(result)


def get_cif_files_generator(input_path):
    if input_path.is_file():
        yield input_path
    else:
        for cif_file in input_path.glob("*.cif"):
            yield cif_file


def process_cif_files_streaming(input_path, output_csv, cif_to_bb_path, decomp_params, num_cores, batch_size, decompose_mode, output_base_dir, verbose, general_timeout, decomp_timeout, obabel_timeout):
    cif_files = list(get_cif_files_generator(input_path))
    total_files = len(cif_files)
    
    if total_files == 0:
        print("No CIF files found")
        return []
    
    if num_cores is None:
        if total_files > 10000:
            num_cores = min(cpu_count() // 2, 8)
        else:
            num_cores = min(cpu_count(), total_files)
    
    print(f"Processing {total_files} files with {num_cores} parallel workers")
    print(f"Per-structure timeout: {general_timeout} seconds")
    
    if decompose_mode and output_base_dir:
        output_base_dir.mkdir(parents=True, exist_ok=True)
    
    valid_count = 0
    processed_count = 0
    timeout_count = 0
    
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cif', 'OxiChecker Validity', 'PATH1_InChI', 'PATH2_InChI', 'PATH3_SMILES', 'Valid_Path'])
        writer.writeheader()
        f.flush()
        
        result_queue = Queue()
        active_processes = {}
        cif_iter = iter(cif_files)
        
        with tqdm(total=total_files, desc="Processing CIF files") as pbar:
            for _ in range(min(num_cores, total_files)):
                try:
                    cif_file = next(cif_iter)
                    p = Process(target=worker_process, 
                              args=(cif_file, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout, result_queue))
                    p.start()
                    active_processes[p] = (cif_file, time.time())
                except StopIteration:
                    break
            
            while active_processes or processed_count < total_files:
                current_time = time.time()
                for p, (cif_file, start_time) in list(active_processes.items()):
                    if current_time - start_time > general_timeout:
                        p.terminate()
                        p.join(timeout=5)
                        if p.is_alive():
                            p.kill()
                            p.join()
                        # Clean CIF name by removing P1_ prefix if present
                        cif_name = Path(cif_file).name
                        if cif_name.startswith('P1_'):
                            cif_name = cif_name[3:]  # Remove 'P1_' prefix
                        timeout_result = {
                            'cif': cif_name,
                            'OxiChecker Validity': 'Processing timeout',
                            'PATH1_InChI': '',
                            'PATH2_InChI': '',
                            'PATH3_SMILES': '',
                            'Valid_Path': ''
                        }
                        writer.writerow(timeout_result)
                        f.flush()
                        
                        processed_count += 1
                        timeout_count += 1
                        pbar.update()
                        
                        del active_processes[p]
                        
                        try:
                            new_cif = next(cif_iter)
                            new_p = Process(target=worker_process,
                                          args=(new_cif, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout, result_queue))
                            new_p.start()
                            active_processes[new_p] = (new_cif, time.time())
                        except StopIteration:
                            pass
                while not result_queue.empty():
                    try:
                        result = result_queue.get_nowait()
                        
                        writer.writerow(result)
                        f.flush()
                        
                        if result['OxiChecker Validity'].startswith('True'):
                            valid_count += 1
                        
                        processed_count += 1
                        pbar.update()
                    except:
                        break
                
                for p in list(active_processes.keys()):
                    if not p.is_alive():
                        p.join(timeout=1)
                        cif_file, start_time = active_processes[p]
                        del active_processes[p]
                        
                        try:
                            new_cif = next(cif_iter)
                            new_p = Process(target=worker_process,
                                          args=(new_cif, cif_to_bb_path, decomp_params, decompose_mode, output_base_dir, verbose, decomp_timeout, obabel_timeout, result_queue))
                            new_p.start()
                            active_processes[new_p] = (new_cif, time.time())
                        except StopIteration:
                            pass
                
                time.sleep(0.1)
            
            for p in active_processes:
                p.terminate()
                p.join(timeout=5)
                if p.is_alive():
                    p.kill()
                    p.join()
    
    print(f"\nResults saved to {output_csv}")
    print(f"Valid: {valid_count}/{total_files}")
    print(f"Timeouts: {timeout_count}/{total_files}")
    
    if decompose_mode and output_base_dir:
        print(f"Decomposition files saved to: {output_base_dir}")
    
    return []


def parse_arguments():
   parser = argparse.ArgumentParser(description='MOF Validity Checker')
   parser.add_argument('--input', required=True, help='CIF file or directory')
   parser.add_argument('--output', default='result.csv', help='Output CSV file')
   parser.add_argument('--version', choices=['v2', 'v3', 'external'], default='v3',
                      help='Choose binary version: v2, v3, or external (default: v3)')
   parser.add_argument('--cif-to-bb-path', help='Path to external cif_to_building_blocks binary (used with --version external)')
   parser.add_argument('--num_cores', type=int, default=None, help='Number of parallel cores')
   parser.add_argument('--batch_size', type=int, default=1000, help='Batch size for processing (deprecated, not used)')
   parser.add_argument('--nosave-decompose-data', action='store_true', help='Do not save decomposition files (delete temporary files)')
   parser.add_argument('--decompose-dir', help='Custom directory for decomposition files (default: ./decompose_<output_name>)')
   parser.add_argument('--verbose', action='store_true', help='Enable verbose logging of decomposition commands')
   
   parser.add_argument('--param1', type=str, default='0.5', help='First decomposition parameter (default: 0.5)')
   parser.add_argument('--param2', type=str, default='5.0', help='Second decomposition parameter (default: 5.0)')
   parser.add_argument('--param3', type=str, default='-0.45', help='Third decomposition parameter (default: -0.45)')
   parser.add_argument('--param4', type=str, default='jmol', help='Fourth decomposition parameter (default: jmol)')
   
   parser.add_argument('--general-timeout', type=int, default=300, help='General timeout per structure in seconds (default: 300)')
   parser.add_argument('--decomposition-timeout', type=int, default=60, help='Decomposition timeout in seconds (default: 60)')
   parser.add_argument('--obabel-timeout', type=int, default=30, help='OpenBabel conversion timeout in seconds (default: 30)')
   
   return parser.parse_args()


def main():
   args = parse_arguments()
   
   global CIF_TO_BB_PATH
   
   if args.version == 'external':
       if not args.cif_to_bb_path:
           print("Error: --cif-to-bb-path required when using --version external")
           sys.exit(1)
       CIF_TO_BB_PATH = Path(args.cif_to_bb_path)
   elif args.version == 'v2':
       CIF_TO_BB_PATH = Path("/opt/libcif/bin/cif_to_building_blocks_v2")
   elif args.version == 'v3':
       CIF_TO_BB_PATH = Path("/opt/libcif/bin/cif_to_building_blocks_v3")
   
   input_path = Path(args.input)
   if not input_path.exists():
       print(f"Error: {input_path} not found")
       sys.exit(1)
   
   output_path = Path(args.output)
   if not output_path.suffix:
       output_path = output_path.with_suffix('.csv')
   
   decomp_params = [args.param1, args.param2, args.param3]
   
   if args.param4:
       decomp_params.append(args.param4)
   
   decompose_mode = not args.nosave_decompose_data
   
   output_base_dir = None
   if decompose_mode:
       if args.decompose_dir:
           output_base_dir = Path(args.decompose_dir)
       else:
           output_name_without_ext = Path(args.output).stem
           output_base_dir = Path(f"./decompose_{output_name_without_ext}")
   
   process_cif_files_streaming(input_path, str(output_path), CIF_TO_BB_PATH, decomp_params, args.num_cores, args.batch_size, decompose_mode, output_base_dir, args.verbose, args.general_timeout, args.decomposition_timeout, args.obabel_timeout)


if __name__ == "__main__":
   main()
