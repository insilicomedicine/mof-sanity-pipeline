# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
from jmolnn_update import JmolNN_update
from ligand import Ligand
import networkx as nx
import logging
logger = logging.getLogger(__name__)


def parse_charges_from_structure_and_graph(s, g):
    N = len(s)
    visited_sites = [False]*N
    metal_idxs = []
    ligands = []
    metal_idxs = [ idx for idx, sp in enumerate(s) if sp.specie.is_metal]
    counter = 1
    if len(metal_idxs) == 0:
        logger.info("It seems to be no metals. Exit OxiChecker check")
        return {"oxichecker_status": "NO METAL", "oxichecker_data" : None}
    try:
        for idx, sp in enumerate(s):
            if not sp.specie.is_metal:
                ligand_data = []
                idx_in_ligands = [x for l in ligands for x in l.idxs]
                Ligand.detect_ligand(idx, g, visited_sites, ligand_data, (0,0,0), idx_in_ligands)
                if len(ligand_data)>0:
                    ligands.append(Ligand(ligand_data, s, counter))
                    counter+=1
    except Exception as e:
        logger.error(e)
        return {"oxichecker_status": "FAILED SEPARATION", "oxichecker_data" : None}
    validity = all([l.validity for l in ligands])
    metal_charges = { idx : 0 for idx in metal_idxs }
    result = []
    if validity:
        try:
            total_charge = 0.0
            for ligand in ligands:
                for grouped_charge in ligand.grouped_charges:
                    total_charge += grouped_charge.detect_metals_charges_and_free_charge(g, metal_idxs, metal_charges)
            result = [(s[x].specie.symbol, y) for x,y in metal_charges.items()]
            result += [("Cell", total_charge)]
        except Exception as e:
            logger.error(e)
            validity = False
    else:
        logger.error("Some of the ligands failed during parsing")
    ligands_data = {}
    for l in ligands:
        mhash = l.graph_hashes['molecule']
        if mhash in ligands_data:
            ligands_data[mhash][0] += 1
        else:
            ligands_data[mhash] = [1, nx.readwrite.json_graph.adjacency_data(l.graph), l.graph_hashes]
    return {"oxichecker_status" : "PASSED",
            "oxichecker_data" : {
                'validity' : validity,
                "ligands_data" : ligands_data,
                "metal_charges" : result
            }
    }
