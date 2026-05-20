# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
from pathlib import Path
import json
import os
import tqdm
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import click
from typing import List, Tuple

equality = lambda x1, x2, tol=0.02: True if abs(x1 - x2)<=tol else False
full_equality = lambda x1, List_x2, tol=0.02: True if any([equality(x1,x2,tol) for x2 in List_x2] ) else False

class Analyzer:

    NORMAL_CHARGES = {'Cu' : [2.0], 'Zn' : [2.0], 'Ni' : [2.0, 3.0],  'Cr' : [3.0, 6.0], 'V' : [3.0, 5.0], 'Au' : [3.0],
                      'Na' : [1.0], 'Cd' : [2.0, 3.0], 'Al' : [3.0], 'Mo' : [4.0, 6.0],  'Fe' : [2.0, 3.0], 'Dy' : [3.0],
                      'K' : [1.0], 'Co' : [2.0, 3.0], 'Ag' : [1.0], 'Pr' : [3.0], 'Sc' : [3.0], 'Sr' : [2.0], 'Ho' : [3.0],
                      'Ca' : [2.0], 'Li' : [1.0], "Zr" : [3.0], 'Gd' : [3.0], 'Pb' : [2.0], 'Ga' : [3.0], 'Rh' : [3.0],
                      'Lu' : [3.0], 'Tm' : [3.0], 'Nb' : [5.0],  'Th' : [4.0], 'Pu' : [4.0], 'Ir' : [3.0, 4.0],
                      'Mn' : [2.0, 3.0], 'Ba' : [2.0], 'Mg' : [2.0], 'Hg' : [2.0], 'Sn' : [0.0, 2.0, 4.0], 'U' : [6.0],
                      "La" : [3.0], 'Yb' : [3.0], 'Nd' : [3.0] , 'Rb' : [1.0], 'Sm' : [3.0], 'Er' : [3.0], 'Eu' : [3.0],
                      'Tb' : [3.0], 'Tl' : [3.0], 'Ce' : [1.0], 'Ti' : [2.0, 4.0], 'Pt' : [2.0, 4.0], 'Ru' : [3.0, 4.0],
                      'Bi' : [3.0], 'In' : [3.0], 'W' : [6.0], 'Cs' : [1.0], 'Pd' : [2.0], 'Hf' : [4.0], "Y" : [3.0]}

    def get_potential_cell_charge(self):
        cell_charges = [0.0]
        for metal in self.metals.keys():
            N = len(cell_charges)
            if len( Analyzer.NORMAL_CHARGES[metal] ) == 1:
                for idx in range(N):
                    cell_charges[idx] += Analyzer.NORMAL_CHARGES[metal][0]*self.metals_counter[metal]
            else:
                cell_charges = cell_charges*len( Analyzer.NORMAL_CHARGES[metal] )
                for idx in range(N):
                    for idx_2 in range(len( Analyzer.NORMAL_CHARGES[metal] )):
                        cell_charges[idx_2*N + idx] += Analyzer.NORMAL_CHARGES[metal][idx_2]*self.metals_counter[metal]
        return cell_charges

    def __init__(self, lines):
        if "Cell" != lines[-1][0]:
            self.validity = False
            return
        self.cell_charge = float(lines[-1][1])
        self.metals = {}
        self.metals_counter = {}
        self.has_unknown_metals = False
        for line in lines[:-1]:
            metal, add_charge = line
            if metal not in Analyzer.NORMAL_CHARGES:
                self.has_unknown_metals = True
            add_charge = float(add_charge)
            adder = self.metals_counter.get(metal, 0) + 1
            charge = self.metals.get(metal, 0.0) + add_charge
            self.metals[metal] = charge
            self.metals_counter[metal] = adder
        for key in self.metals:
            self.metals[key] /= self.metals_counter[key]
        self.validity = True

    def check(self):
        if not self.validity:
            self.good = 3
            return self.good
        if self.has_unknown_metals:
            self.good = 2
            return self.good

        cell_charge = 0.0
        counter = 0
        for key in self.metals:
            cell_charge += self.metals[key]*self.metals_counter[key]
            counter += self.metals_counter[key]
        cell_charge -= self.cell_charge
        res = self.get_potential_cell_charge()
        if full_equality(cell_charge, res):
            self.good = 1
            return self.good
        else:
            self.good = 0
            return self.good


def read_json(json_file):
    with open(json_file) as f:
        data = json.load(f)
    return data

def single_json_parser(result_dir_path, cif_stem):
    try:
        # New structure: JSON files are in result_dir_path with name cif_stem.json
        json_file_path = os.path.join(result_dir_path, cif_stem + '.json')
        
        if os.path.exists(json_file_path):
            data = read_json(json_file_path)
        else:
            # Try the old structure with suffix files
            data = {}
            json_file_no_ext = os.path.join(result_dir_path, cif_stem)
            for suffix in ['base', 'libcif', 'oxichecker', 'framechecker', 'platon']:
                suffix_file = json_file_no_ext + "__" + suffix + '.json'
                if os.path.exists(suffix_file):
                    data.update(read_json(suffix_file))
        if len(data)==0:
            return (None, None, "EMPTY_DATA")
        tmp_dict = {}
        for key in  ['content_hash', 'basic_validity',
                    'formula', 'reduced_formula', 'density', 'volume',
                    'group_str', 'group_id', 'structure_hash_strict',
                    'structure_hash',
                    'is_graph_constructed',
                    'graph_dim', 'platon_validity',
                    'libcif_validity', 'framecheker_validity']:
            if key in data:
                tmp_dict[key] = data[key]
            if key == 'basic_validity':
                if not isinstance(data[key], bool):
                    print(data)
        if 'platon_validity' in tmp_dict and tmp_dict['platon_validity']:
            platon_bad_lines = len([x for x in data['platon_data'] if "Unusual sp" in x])
            if platon_bad_lines>0:
                tmp_dict['platon_condition_is_good'] = False
            else:
                tmp_dict['platon_condition_is_good'] = True
        else:
            tmp_dict['platon_condition_is_good'] = None

        if 'libcif_validity' in tmp_dict and tmp_dict['libcif_validity']:
            libcif_out = data.get('libcif_out', '').strip()
            if libcif_out:
                res = libcif_out.split()[0]
                tmp_dict['libcif_result'] = True if res=='OK' else False
                if tmp_dict['libcif_result']==False:
                    libcif_error = data.get('libcif_error', '')
                    if libcif_error and 'graph' in libcif_error.lower():
                        tmp_dict['libcif_result'] = True
            else:
                tmp_dict['libcif_result'] = False

            libcif_error = data.get('libcif_error', '') or ''
            libcif_warning_lines = []
            for line in libcif_error.splitlines():
                normalized = line.strip()
                if not normalized:
                    continue
                upper_line = normalized.upper()
                if 'BAD' in upper_line or 'WARNING' in upper_line:
                    libcif_warning_lines.append(normalized)

            if libcif_warning_lines:
                tmp_dict['libcif_warning'] = ' | '.join(libcif_warning_lines)
                if tmp_dict.get('libcif_result') is True:
                    tmp_dict['libcif_result'] = 'Inspect'
            else:
                tmp_dict['libcif_warning'] = None

        else:
            tmp_dict['libcif_validity'] = False


        if 'framecheker_validity' in tmp_dict and tmp_dict['framecheker_validity']:
            tmp_dict.update(data['framecheker_results'])

        else:
            tmp_dict['framecheker_validity'] = False

        return (cif_stem, tmp_dict, "SUCCESS")
    except Exception as e:
        return (None, None, f"EXCEPTION_{type(e).__name__}_{str(e)[:50]}")


def find_valid_cifs(folder: str) -> List[Tuple[str, str]]:
    """
    Collects all result directories from the new output structure.
    Looks for individual CIF result directories in results/ subdirectory.

    Args:
        folder (str): Path to the output directory containing results/ subdirectory.

    Returns:
        List[Tuple[str, str]]: List of (results directory path, CIF stem) tuples.
    """
    folder_path = Path(folder)
    results_dir = folder_path / "results"
    
    if not results_dir.exists():
        # Fallback: try old structure - look for CIF files directly
        return [
            (str(x.parent), x.stem)
            for x in folder_path.glob('*.cif')
            if "OBABEL_" not in x.name and "P1_" not in x.name and "SYMM_" not in x.name
        ]
    
    # New structure: look for result directories
    result_dirs = []
    for subdir in results_dir.iterdir():
        if subdir.is_dir():
            # Check if this directory has JSON files
            json_files = list(subdir.glob('*.json'))
            if json_files:
                result_dirs.append((str(subdir), subdir.name))
    
    return result_dirs


def format_cif_name(path_str):
    """
    Format the CIF file path for the 'cif' column.
    Removes './' prefix and adds '.cif' extension if needed.
    """
    if path_str.startswith('./'):
        path_str = path_str[2:]
    if not path_str.endswith('.cif'):
        path_str = path_str + '.cif'
    return path_str


@click.command()
@click.option('--folder', '-f', type=click.Path(exists=True, file_okay=False), default='.',
              help='Output folder containing results/ subdirectory.', show_default=True)
@click.option('--n-jobs', '-j', type=int, default=50, help='Number of parallel jobs.', show_default=True)
@click.option('--output', '-o', type=click.Path(dir_okay=False), default='sanity_results.csv',
              help='CSV file to save results.', show_default=True)
@click.option('--oxichecker-path', '-m', type=click.Path(exists=True, dir_okay=False), default=None,
              help='Path to OxiChecker CSV file to merge with results.')
def main(folder: str, n_jobs: int, output: str, oxichecker_path: str):
    """
    CLI tool to process sanity runner results and save to a CSV file.
    Expects the new output structure with results/ subdirectory containing individual CIF result folders.
    """
    cifs = find_valid_cifs(folder)
    click.echo(f"Found {len(cifs)} result directories in '{folder}'.")

    results = Parallel(n_jobs=n_jobs, verbose=7)(delayed(single_json_parser)(result_dir_path, cif_stem) for result_dir_path, cif_stem in cifs)
    
    # Debug output: analyze parsing results
    successful_results = [(i, j) for i, j, status in results if i is not None]
    failed_results = [(status,) for i, j, status in results if i is None]
    
    full_data = {i:j for i,j in successful_results}
    
    # Count error types
    error_counts = {}
    for (error_type,) in failed_results:
        error_counts[error_type] = error_counts.get(error_type, 0) + 1
    
    click.echo(f"\n=== JSON PARSING RESULTS ===")
    click.echo(f"Successfully parsed: {len(successful_results)}/{len(results)} structures")
    click.echo(f"Failed to parse: {len(failed_results)}/{len(results)} structures")
    if error_counts:
        click.echo("Error breakdown:")
        for error_type, count in sorted(error_counts.items()):
            click.echo(f"  {error_type}: {count}")
    click.echo("=" * 29)

    df = pd.DataFrame(full_data).T
    
    df.index = df.index.map(format_cif_name)
    df.index.name = 'cif'
    df = df.reset_index()

    df = df.rename(columns = {'basic_validity' : "Basic Validity", "libcif_result" : "LibCif Validity",
"platon_condition_is_good" : "PLATON Validity", 'libcif_warning': 'LibCif_Warning'})

    for validity_col in ["Basic Validity", "PLATON Validity"]:
        if validity_col in df.columns:
            df[validity_col] = df[validity_col].astype(bool).fillna("Check not performed")
        else:
            df[validity_col] = "Check not performed"

    if "LibCif Validity" in df.columns:
        def _normalize_libcif_validity(value):
            if pd.isna(value):
                return "Check not performed"  
            if isinstance(value, str):
                lower_value = value.strip().lower()  
                if lower_value in {"true", "ok"}:  
                    return True
                if lower_value in {"false", "fail"}: 
                    return False 
                if lower_value == "inspect":
                    return True
                return value 
            return bool(value)
        df["LibCif Validity"] = df["LibCif Validity"].apply(_normalize_libcif_validity)
    else:
        df["LibCif Validity"] = "Check not performed"

    if "LibCif_Warning" not in df.columns:
        df["LibCif_Warning"] = np.nan

    checker_keys = ['HAS_ATOMIC_OVERLAPS', 'HAS_OVERCOORDINATED_C', 'HAS_OVERCOORDINATED_H',
                    'HAS_OVERCOORDINATED_N', 'HAS_LONE_MOLECULE', 'HAS_BAD_RARE_EARTH',
                    'HAS_BAD_ALKALI_ALKALINE', 'HAS_BAD_TERMINAL_OXO']

    for key in checker_keys:
        renamer = {x : x.lower().replace("_", " ") for x in checker_keys}
        df = df.rename(columns = renamer)
        df[renamer[key]] = df[renamer[key]].astype(bool).fillna("Check not performed")
    
    # Save intermediate CSV before OxiChecker merge
    intermediate_output = output.replace('.csv', '_pre_oxichecker.csv')
    df.to_csv(intermediate_output, index=False)
    click.echo(f"Intermediate results (before OxiChecker merge) saved to {intermediate_output}")

    if oxichecker_path:
        click.echo(f"Loading OxiChecker data from '{oxichecker_path}'...")
        try:
            oxichecker_df = pd.read_csv(oxichecker_path)
            click.echo(f"OxiChecker CSV contains {len(oxichecker_df)} rows and columns: {list(oxichecker_df.columns)}")
            
            # Normalize CIF column name - handle different possible column names
            cif_column = None
            for col in ['cif', 'CIF', 'filename', 'file', 'structure']:
                if col in oxichecker_df.columns:
                    cif_column = col
                    break
            
            if cif_column is None:
                click.echo("Warning: No CIF filename column found in OxiChecker file. Using first column.")
                cif_column = oxichecker_df.columns[0]
            
            # Rename to standard 'cif' column name
            if cif_column != 'cif':
                oxichecker_df = oxichecker_df.rename(columns={cif_column: 'cif'})
            
            # Normalize CIF filenames for matching
            oxichecker_df['cif'] = oxichecker_df['cif'].apply(format_cif_name)
            
            # Check for validity column
            validity_column = None
            for col in ['OxiChecker Validity', 'validity', 'valid', 'is_valid', 'decomposed']:
                if col in oxichecker_df.columns:
                    validity_column = col
                    break
            
            if validity_column and validity_column != 'OxiChecker Validity':
                oxichecker_df = oxichecker_df.rename(columns={validity_column: 'OxiChecker Validity'})
            
            # Merge OxiChecker data with sanity results
            original_count = len(df)
            df = pd.merge(df, oxichecker_df, on='cif', how='left')
            matched_count = df['OxiChecker Validity'].notna().sum()
            
            click.echo(f"Successfully merged OxiChecker data: {matched_count}/{original_count} structures matched")
            
            # If no validity column found, add a default one
            if 'OxiChecker Validity' not in df.columns:
                click.echo("Warning: No validity column found in OxiChecker file. Adding 'OxiChecker Validity' as 'Processed'")
                df['OxiChecker Validity'] = 'Processed'
                
        except Exception as e:
            click.echo(f"Error loading OxiChecker data: {e}")
            df['OxiChecker Validity'] = 'Error loading OxiChecker data'
    else:
        df['OxiChecker Validity'] = 'Not processed'

    df = df [ ['cif'] + [x for x in ['content_hash', 'formula', 'reduced_formula',
        'density', 'volume', 'group_str', 'group_id', 'structure_hash_strict',
        'structure_hash', 'is_graph_constructed', 'graph_dim',  'HAS_OMS',
        'DECORATED_GRAPH_HASH', 'UNDECORATED_GRAPH_HASH',
        'DECORATED_SCAFFOLD_HASH', 'UNDECORATED_SCAFFOLD_HASH',
            'Basic Validity', 'LibCif Validity', 'LibCif_Warning', 'PLATON Validity', 'OxiChecker Validity',
        'has atomic overlaps', 'has overcoordinated c', 'has overcoordinated h',
        'has overcoordinated n', 'has lone molecule', 'has bad rare earth',
        'has bad alkali alkaline', 'has bad terminal oxo'] if x in df.columns ] ]

    main_validity_cols = ["Basic Validity", "LibCif Validity", "PLATON Validity"]
    for col in main_validity_cols:
        if col in df.columns and df[col].dtype == object:
            df.loc[df[col] == "Check not performed", col] = False
    
    if 'OxiChecker Validity' in df.columns:
        oxichecker_good = df['OxiChecker Validity'].astype(str).str.contains('True', na=False)
    else:
        oxichecker_good = pd.Series([False] * len(df), index=df.index)

    ok_cond = np.logical_and.reduce(
        [df[key]!=True for key in [x.lower().replace("_", " ") for x in checker_keys] if key in df.columns] + \
        [df[key]!=False for key in main_validity_cols if key in df.columns] + \
        [np.logical_or(oxichecker_good, df['OxiChecker Validity'] == 'Not processed')]
    )
    
    df.loc[:, 'Sanity'] = ok_cond
    
    # Reorder columns to put validation results first
    validation_cols = ['cif', 'Basic Validity', 'LibCif Validity', 'LibCif_Warning', 'PLATON Validity', 'OxiChecker Validity', 'Sanity']
    existing_validation_cols = [col for col in validation_cols if col in df.columns]
    other_cols = [col for col in df.columns if col not in existing_validation_cols]
    df = df[existing_validation_cols + other_cols]
    
    df.to_csv(output, index=False)

    # Print summary statistics
    total_structures = len(df)
    basic_valid = (df['Basic Validity'] == True).sum() if 'Basic Validity' in df.columns else 0
    libcif_valid = (df['LibCif Validity'] == True).sum() if 'LibCif Validity' in df.columns else 0
    platon_valid = (df['PLATON Validity'] == True).sum() if 'PLATON Validity' in df.columns else 0
    
    # OxiChecker validity counting
    oxichecker_valid = 0
    oxichecker_processed = 0
    if 'OxiChecker Validity' in df.columns:
        oxichecker_processed = (df['OxiChecker Validity'] != 'Not processed').sum()
        oxichecker_valid = df['OxiChecker Validity'].astype(str).str.contains('True', na=False).sum()
    
    sanity_passed = (df['Sanity'] == True).sum() if 'Sanity' in df.columns else 0
    
    click.echo(f"\nResults saved to {output}")
    click.echo(f"\n=== SUMMARY STATISTICS ===")
    click.echo(f"Total structures: {total_structures}")
    click.echo(f"Basic validity: {basic_valid}/{total_structures} ({100*basic_valid/total_structures:.1f}%)")
    click.echo(f"LibCif validity: {libcif_valid}/{total_structures} ({100*libcif_valid/total_structures:.1f}%)")
    click.echo(f"PLATON validity: {platon_valid}/{total_structures} ({100*platon_valid/total_structures:.1f}%)")
    click.echo(f"OxiChecker processed: {oxichecker_processed}/{total_structures} ({100*oxichecker_processed/total_structures:.1f}%)")
    click.echo(f"OxiChecker valid: {oxichecker_valid}/{oxichecker_processed if oxichecker_processed > 0 else 1} ({100*oxichecker_valid/(oxichecker_processed if oxichecker_processed > 0 else 1):.1f}%)")
    click.echo(f"Overall sanity: {sanity_passed}/{total_structures} ({100*sanity_passed/total_structures:.1f}%)")
    click.echo("=" * 27)


if __name__ == '__main__':
    main()
