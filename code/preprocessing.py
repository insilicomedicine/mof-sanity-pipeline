#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
'''Preprocessing of the initial cif (or other files) to the unified pymatgen json format.
This module is separated in order to hide the realiztion of the cif 
'''
import hashlib
import logging
import os
import numpy as np
from typing import Optional
from jmolnn_update import JmolNN_update
from openbabel import openbabel
from pymatgen.io.cif import CifParser
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.analysis.structure_matcher import ElementComparator
from pymatgen.core import Structure
import warnings
from openbabel import pybel
import signal
from contextlib import contextmanager

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
pybel.ob.obErrorLog.StopLogging()

SYMPREC = 0.01
ANGLE_TOL = 5
JMOLNN_TOLERANCE = 0.45
JMOLNN_MIN_DIST = 0.001
STRUCTURE_PARSE_TIMEOUT = 60  
GRAPH_BUILD_TIMEOUT = 60     

PROPLIST = ['content_hash',
'basic_validity',
'matgen_json', 
'formula', 
'reduced_formula',
'density', 
'volume', 
'group_str', 
'group_id', 
'structure_hash_strict', 
'structure_hash',
'sym_filepath']

GRAPH_PROPS = ['is_graph_constructed', 
               'graph_json', 'graph_dim']

logger = logging.getLogger(__name__)

@contextmanager
def timeout(seconds):
    def timeout_handler(signum, frame):
        raise TimeoutError(f"Operation timed out after {seconds} seconds")
    
    old_handler = signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

class FileStructure:
    _extension : Optional[str] = None
    _content : Optional[str] = None
    
    content_hash : Optional[str] = None
    basic_validity : bool = False
    formula : Optional[str] = None
    volume : Optional[float] = None
    density : Optional[float] = None
    group_str : Optional[str] = None
    group_id : Optional[int] = None
    matgen_structure = None
    p1_filepath : Optional[str] = None
    matgen_json : Optional[str] = None
    reduced_formula : Optional[dict] = None
    structure_hash : Optional[str] = None
    structure_hash_strict : Optional[str] = None
    is_graph_constructed : Optional[bool] = False
    graph = None
    graph_json : Optional[dict] = None
    graph_dim : Optional[dict] = None
    sym_filepath : Optional[str] = None


    @property
    def content(self):
        if self._content is None:
            logger.info(f"Reading content")
            try:
                with open(self.filepath) as f:
                    self._content = f.read().strip().replace("'",'"')
            except:
                logger.error(f"Could not open file or read content")
                raise Exception
            logger.info(f"{len(self._content)} symbols read")
        return self._content
    
    @property
    def file_content_hash(self):
        if self.content_hash is None:
            self.content_hash =  hashlib.md5(self.content.encode()).hexdigest()
            logger.info(f"Calculated striped file hash = {self.content_hash}")
        return self.content_hash
    
    @property
    def extension(self):
        if self._extension is None:
            self._extension = str(self.filepath).split(".")[-1].lower()
            logger.info(f"Working with {self.extension.upper()} filetype")
        return self._extension

    def obabel_convertion(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(self.extension, "cif")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f"{self.filepath}") 
        
        # Use babel_cifs directory if available, otherwise use original dirname
        babel_dir = self.cif_dirs.get('babel_cifs', self.dirname)
        self.filepath = os.path.join(str(babel_dir), 
                                     f"OBABEL_{self.basename.replace(f'.{self.extension}', '.cif')}")
        self._extension = 'cif'
        obConversion.WriteFile(mol, self.filepath)
        logger.info(f"Obabel conversion performed")
        logger.info(f"Updated filepath to {self.filepath}")

    @staticmethod
    def make_structure_hash(matgen_structure, regime = 'mild'):
        if regime == 'mild':
            atom_array = np.array([[a,x,y,z] for a,(x,y,z) in zip(matgen_structure.atomic_numbers, 
                                      (np.round(matgen_structure.frac_coords, 2)*10**2).astype(int))])
            position_string = " ".join(map(str, [atom_array.shape[0]] + \
                        list( (np.round(matgen_structure.lattice.abc, 1)*10**1).astype(int)) + \
                        list((np.round(matgen_structure.lattice.angles, 0)).astype(int)) + \
                        list(atom_array[np.lexsort(np.transpose(atom_array)[::-1])].ravel())))
        else:
            atom_array = np.array([[a,x,y,z] for a,(x,y,z) in zip(matgen_structure.atomic_numbers, 
                                       (np.round(matgen_structure.frac_coords, 3)*10**3).astype(int))])
            position_string = " ".join(map(str, [atom_array.shape[0]] + \
                                    list((np.round(matgen_structure.lattice.abc, 2)*10**3).astype(int)) + \
                                    list((np.round(matgen_structure.lattice.angles, 1)*10**1).astype(int)) + \
                                    list(atom_array[np.lexsort(np.transpose(atom_array)[::-1])].ravel())))
        
        structure_hash = hashlib.md5(position_string.encode("utf-8")).hexdigest()
        return structure_hash
    
    def prepare_clean_struture(self,
                               matgen_symprec : float = SYMPREC, 
                               matgen_angle_tolerance : float = ANGLE_TOL):
        warnings.filterwarnings("ignore")
        try:
            file_size_mb = os.path.getsize(self.filepath) / (1024 * 1024)
            if file_size_mb > 50:  
                logger.error(f"File too large: {file_size_mb:.2f} MB, skipping")
                self.basic_validity = False
                return
            
            logger.info(f"Starting structure parsing with {STRUCTURE_PARSE_TIMEOUT}s timeout")
            with timeout(STRUCTURE_PARSE_TIMEOUT):
                self.matgen_structure = CifParser.from_str(self.content).parse_structures(primitive = True)[0]
            
            logger.info(f"Primitive structure calculated")
            
            num_atoms = len(self.matgen_structure)
            if num_atoms > 10000:  
                logger.error(f"Structure too large: {num_atoms} atoms, skipping")
                self.basic_validity = False
                self.matgen_structure = None  
                return
            
            self.formula = self.matgen_structure.formula
            self.reduced_formula = ElementComparator().get_hash(self.matgen_structure.composition).as_dict()
            self.density = float(self.matgen_structure.density)
            self.volume = self.matgen_structure.volume
            logger.info(f"Formula: {self.formula}, Volume: {self.volume}, Density: {self.density}")
            self.group_str, self.group_id = \
                    self.matgen_structure.get_space_group_info(matgen_symprec, matgen_angle_tolerance)
            self.structure_hash_strict = FileStructure.make_structure_hash(self.matgen_structure, 'strict')
            self.structure_hash = FileStructure.make_structure_hash(self.matgen_structure, 'mild')
            logger.info(f"Nonstrict 3D structure hash: {self.structure_hash}")
            
            # Use p1_cifs directory if available, otherwise use original dirname
            p1_dir = self.cif_dirs.get('p1_cifs', self.dirname)
            self.p1_filepath = os.path.join(str(p1_dir), 
                                    f"P1_{self.basename}")
            _ = self.matgen_structure.to(self.p1_filepath, fmt='cif', symprec=None)
            logger.info(f"P1 primitive structure written")
            logger.info(f"Group is: {self.group_id}, {self.group_str}")
            if self.group_id != 1:
                # Use symm_cifs directory if available, otherwise use original dirname
                symm_dir = self.cif_dirs.get('symm_cifs', self.dirname)
                self.sym_filepath = os.path.join(str(symm_dir), 
                                    f"SYMM_{self.basename}")
                _ = self.matgen_structure.to(self.sym_filepath, 
                                            fmt='cif', 
                                            symprec = matgen_symprec,
                                            angle_tolerance = matgen_angle_tolerance)
                logger.info(f"SYMM primitive structure written")
            else:
                self.sym_filepath = None
            self.matgen_json = self.matgen_structure.to_json()
            self.basic_validity = True
            logger.info(f"All based properties have been calculated. Structure is valid")
        except TimeoutError as e:
            logger.error(f"Structure parsing timeout: {e}")
            self.basic_validity = False
            self.matgen_structure = None  
        except Exception as e:
            logger.error(f"Matgen Structure Creation Error: {e}")
            self.basic_validity = False

    def prepare_graph(self):
        try:
            logger.info(f"Starting graph construction with {GRAPH_BUILD_TIMEOUT}s timeout")
            with timeout(GRAPH_BUILD_TIMEOUT):
                self.graph = JmolNN_update(self.matgen_structure,
                                           tol = JMOLNN_TOLERANCE, 
                                           min_bond_distance = JMOLNN_MIN_DIST).\
                                            get_bonded_structure(self.matgen_structure)
            
            logger.info(f"Structure Graph Constructed")
            self.graph_json = self.graph.to_json()
            self.graph_dim = int(get_dimensionality_larsen(self.graph))
            self.is_graph_constructed = True
        except TimeoutError as e:
            logger.error(f"Graph construction timeout: {e}")
            self.is_graph_constructed = False
            self.graph = None 
        except Exception as e:
            logger.error(f"Graph Creation Error: {e}")
            self.is_graph_constructed = False
    
    def __init__(self, filepath, obabel_run, cif_dirs=None):
        self.filepath = str(filepath)
        self.basename = os.path.basename(self.filepath)
        self.dirname = os.path.dirname(self.filepath)
        self.obabel_run = obabel_run
        self.cif_dirs = cif_dirs or {}
        if obabel_run:
            self.obabel_convertion()
        assert self.extension == 'cif'
        self.basename = '.'.join(self.basename.split('.')[:-1]) + '.cif'
        self.content
        self.file_content_hash
        self.prepare_clean_struture()
        self.prepare_graph()
    
    def property_dict(self):
        if self.basic_validity:
            data = {x:self.__dict__[x] for x in PROPLIST}
            if self.is_graph_constructed:
                data.update( {x:self.__dict__[x] for x in GRAPH_PROPS} )
            else:
                data.update( {"is_graph_constructed" : False } )
            return data
        else:
            return {x:self.__dict__[x] for x in ['content_hash', 'basic_validity']}