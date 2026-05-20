#!/usr/bin/env python
# coding: utf-8
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.

import subprocess
import os
from joblib import Parallel, delayed
from typing import List, Literal, Generator
import pandas as pd
import warnings
import argparse
import pathlib
from pathlib import Path
import time
import json

R_MIN_DEFAULT = 0.5
R_MAX_DEFAULT = 5.0
STRATEGY_DEFAULT = 'cordero'
STRATEGY_CHOISES = ['cordero', 'pyykko']
TOLERANCE_DEFAULT = 1.15
BIN_CHOISES = ['cif_check', 'cif_ok']
BIN_DEFAULT = 'cif_check'

def get_valid_kwargs(func: callable, args_dict: dict) -> dict:
    """
    Filters and returns valid keyword arguments for a given function.

    Parameters:
        func (callable): The function to check for valid keyword arguments.
        args_dict (dict): The dictionary of arguments to filter.

    Returns:
        dict: A dictionary containing only the valid keyword arguments for the function.
    """
    # Get the function's argument names and count of arguments
    valid_args = func.__code__.co_varnames[:func.__code__.co_argcount]
    kwargs_len = len(func.__defaults__) if func.__defaults__ else 0  # Number of keyword arguments
    valid_kwargs = valid_args[-kwargs_len:]  # Keyword arguments are the last ones

    # Filter and return the valid keyword arguments
    return {key: value for key, value in args_dict.items() if key in valid_kwargs}

def chunks(xs : list, n : int) -> Generator:
    n = max(1, n)
    return (xs[i:i+n] for i in range(0, len(xs), n))

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')

def get_cifs_for_folder(folder : str, recursive : bool = False) -> List[pathlib.Path]:
    if recursive:
        return [file for file in list(Path(folder).rglob("*.[cC][iI][fF]")) if file.is_file()]
    else:
        return [file for file in list(Path(folder).glob("*.[cC][iI][fF]")) if file.is_file()]

def parse_cif_input(cifs_path : List[str], recursive : bool = False) -> List[pathlib.PosixPath]:
    parsed_paths = []
    for str_path in cifs_path:
        if os.path.isfile(str_path):
            if str_path.lower().endswith('.cif'):
                parsed_paths.append(pathlib.Path(str_path))
            else:
                if not str_path.lower().endswith('.cif') and not str_path.lower().endswith('.txt') :
                    warnings.warn(f"File {str_path} has no '.cif' or '.txt' extension. It would be proccessed as txt, but please check the input")
                with open(str_path) as f:
                    parsed_paths += [pathlib.Path(filename) for filename in f.read().split() if filename.lower().endswith('.cif')]
        if os.path.isdir(str_path):
            parsed_paths += get_cifs_for_folder(str_path, recursive=recursive)
    return parsed_paths

class Executable:
    def __init__(self, libcif_folder : str,
                 executable : Literal[*BIN_CHOISES] = BIN_DEFAULT,
                 r_min : float = R_MIN_DEFAULT,
                 r_max : float = R_MAX_DEFAULT,
                 tolerance : float = TOLERANCE_DEFAULT,
                 strategy : Literal[*STRATEGY_CHOISES] = STRATEGY_DEFAULT):
        if not os.path.exists(libcif_folder):
            raise Exception(f"There is no folder: {libcif_folder}")
        if not os.path.isdir(libcif_folder):
            raise Exception(f"{libcif_folder} is not a folder")
        if executable not in BIN_CHOISES:
            raise Exception(f"Sorry. Now only {BIN_CHOISES} executables are parsed")
        self.executable = executable
        self.executable_path = os.path.join(libcif_folder,executable)
        if not os.path.exists(self.executable_path):
            raise Exception(f"There is no such file: {self.executable_pathy}")
        if not os.path.isfile(self.executable_path):
            raise Exception(f"{self.executable_path} is not a file")
        if not os.access(self.executable_path, os.X_OK):
            raise Exception(f"{self.executable_path} is not executable.\nRun 'chmod +x {self.executable_path}' before running this script.")
        if tolerance<1.0001:
            warnings.warn(f"Tolerance value {tolerance:.5f} is too small. Recommended values are between 1.1 and 1.3")
        if tolerance>1.5:
             warnings.warn(f"Tolerance value {tolerance:.5f} is too big. Recommended values are between 1.1 and 1.3")
        self.tolerance = tolerance
        if r_min<0.45:
            warnings.warn(f"R_min value {r_min:.5f} is too small. Recommended value is 0.5")
        if r_min>0.6:
            warnings.warn(f"R_min value {r_min:.5f} is too big. Recommended value is 0.5")
        self.r_min = r_min
        self.r_max = r_max
        self.strategy = strategy
        
    @property
    def list_args(self):
        return self.executable_path, [self.r_min, self.r_max, self.tolerance, self.strategy]

    @staticmethod
    def cif_sanity_out_preprocess(outstring : str) -> dict:
        status, formula, *cell_params = outstring.strip().split()
        Volume_A3, density_g_cm3, a, b, c, alpha, beta, gamma = map(float, cell_params)
        formula = formula.split('_')
        formula = {element : amount for element, amount in zip(formula[::2], formula[1::2])}
        return {"status" : status, "formula" : formula,
                'Volume [A^3]' : Volume_A3, 'density [g/cm^3]' : density_g_cm3,
                'a' : a, 'b' : b, 'c' : c, 'alpha' : alpha, 'beta' : beta, 'gamma': gamma}

    def single_file_run(self, file : pathlib.Path, recursive : bool = False, keep_folder_names : bool = False) -> dict:
        out_dict = {'name': file.name}
        if keep_folder_names:
            out_dict.update({'folder': str(file.parent) if str(file.parent)!='.' else None})
        cpp_binary, arguments = self.list_args
        try:
            runner = subprocess.Popen([cpp_binary, str(file), *[str(arg) for arg in arguments]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = map(lambda x: x.decode('utf-8'), runner.communicate())
            if self.executable == 'cif_check':
                out_dict.update(self.cif_sanity_out_preprocess(stdout))
            else:
                out_dict.update({'status' : stdout.strip()})
            out_dict.update({"error" : stderr.strip()})
        except:
           pass       
        return out_dict


    def run_all(self, paths : List[str], out_file : str,
                recursive : bool = False, keep_folder_names : bool = False,
                n_jobs : int = -1, verbose : int = 5, update_out : int = 0) -> pd.DataFrame:
        results_array = []
        paths = parse_cif_input(paths, recursive=recursive)
        if os.path.exists(out_file):
            new_out_file = '_'.join([str(x) for x in list(time.localtime())[:6]]) + "_" + out_file
            warnings.warn(f"The file {out_file} already exists. Use time mark. The new out file is {new_out_file}")
            out_file = new_out_file
        suffix = out_file.split('.')[-1]
        if (suffix.lower() != 'csv') and (suffix.lower() != 'parquet') and (suffix.lower() != 'json'):
            warnings.warn(f"CSV / PARQUET / JSON are allowed. Will write json")
            suffix = 'json'
        write_out = lambda df : getattr(df, f'to_{suffix}')
        if n_jobs==1 or len(paths)==1:
            print("Serial run")
            for counter, path in enumerate(paths, 1):
                results_array.append(self.single_file_run(path, recursive=recursive, keep_folder_names=keep_folder_names))
                if update_out != 0 and counter%update_out==0:
                    write_out(pd.DataFrame(results_array))(out_file)
            df = pd.DataFrame(results_array)
            write_out(df)(out_file)
        else:
            print('Parallel run')
            if update_out == 0:
                results_array = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(self.single_file_run)(path, recursive=recursive, 
                                                                                                        keep_folder_names=keep_folder_names) for path in paths)
                df = pd.DataFrame(results_array)
                write_out(df)(out_file)
            else:
                if update_out<500:
                    warnings.warn(f"Please increase the --update-out parameter. Now it is too small. Recommended value is more than 1000")
                if update_out<n_jobs:
                    update_out = n_jobs
                    warnings.warn(f"Have changed the --update-out to be equal to --n-jobs {n_jobs}")
                for subarray in chunks(paths, update_out):
                    results_array += Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(self.single_file_run)(path, recursive=recursive, 
                                                                                                        keep_folder_names=keep_folder_names) for path in subarray)
                    write_out(pd.DataFrame(results_array))(out_file)
                df = pd.DataFrame(results_array)
                write_out(df)(out_file)
        return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(f'''Run libcif executable in parallel''', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=120))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)
    required.add_argument('--cifs-path', '-i', required=True, nargs='+', type=str, 
                          help=f"Path to cif file(s) or(and) folder(s) containing cif files. TXT files with cif paths could be added")
    required.add_argument('--out-file', '-o', required=True, nargs=1, type=str,
                          help=f"Output filename to store the result. Must be csv or parquet")
    required.add_argument('--libcif-folder', '-f', required=True, nargs=1, type=str, 
                          help=f"Libcif folder with executables")
    optional.add_argument('--executable', '-x', required=False, nargs='?', type=str, 
                          choices=[*BIN_CHOISES], default=BIN_DEFAULT, help=f"Libcif executable to run. Default to {BIN_DEFAULT}")
    optional.add_argument('--r-min', required=False, nargs='?', type=float, default = R_MIN_DEFAULT, 
                          help=f"The intersection distance to check. Default to {R_MIN_DEFAULT}")
    optional.add_argument('--r-max', required=False, nargs='?', type=float, default = R_MAX_DEFAULT,
                          help=f"The max range distance. Default to {R_MAX_DEFAULT}. Legacy option")
    optional.add_argument('--tolerance', required=False, nargs='?', type=float, default = TOLERANCE_DEFAULT,
                          help=f"The coefficient where there is bond IF r < TOLERANCE*(cov_rad_1 + cov_rad_2). Default to {TOLERANCE_DEFAULT}")
    optional.add_argument('--strategy', required=False, nargs='?', type=str, choices=[*STRATEGY_CHOISES], 
                          default=STRATEGY_DEFAULT, help=f"The strategy to calculate atomic radii. Default to {STRATEGY_DEFAULT}")
    optional.add_argument('--recursive', required=False, type=str_to_bool, nargs='?', const=True, default=False,
                          help=f"Flag to recursively run the script inside the CIFS_PATH folders. Default to false")
    optional.add_argument('--keep-folder-names', required=False, type=str_to_bool, nargs='?', const=True, default=False,
                          help=f"Flag to keep CIFS_PATH folder names together with cif names. Default to false")
    optional.add_argument('--n-jobs', required=False, nargs='?', type=int, default=1, 
                          help=f"Number of CPUs to run. '-1' means run using all cores. Default to 1")
    optional.add_argument('--verbose', required=False, nargs='?', type=int, default=5, 
                          help=f"Verbosity level of joblib [0 to 10]. Default to 5")
    optional.add_argument('--update-out', required=False, nargs='?', type=int, default=0, 
                          help=f"Number of iterations after which one should update the out file. Default to 0 which means no updates till the end")
    args = vars(parser.parse_args())

    executable = Executable(args['libcif_folder'][0], **get_valid_kwargs(Executable.__init__, args))
    executable.run_all(args['cifs_path'], args['out_file'][0], **get_valid_kwargs(Executable.run_all, args))
