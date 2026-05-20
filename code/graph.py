# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
from jmolnn_update import JmolNN_update

JMOLNN_TOLERANCE = 0.45
JMOLNN_MIN_DIST = 0.001

def calculate_graph(structure):
    try:
        graph = JmolNN_update(structure, tol = JMOLNN_TOLERANCE, 
                              min_bond_distance = JMOLNN_MIN_DIST).get_bonded_structure(structure)
        return graph
    except:
        return None