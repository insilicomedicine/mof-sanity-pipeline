# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Module for radii lookups."""
import warnings

import numpy as np

from .definitions import COVALENT_RADII, VDW_RADII

_COVALENT_RADII_MEDIAN = np.median(list(COVALENT_RADII.values()))
_VDW_RADII_MEDIAN = np.median(list(VDW_RADII.values()))


def _get_covalent_radius(element):
    try:
        radius = COVALENT_RADII[element]

    except KeyError:
        radius = _COVALENT_RADII_MEDIAN
        warnings.warn(f"Covalent radius for {element} unknown. Using median {radius:.2f}.")
    return radius


def _get_vdw_radius(element):
    try:
        radius = VDW_RADII[element]

    except KeyError:
        radius = _VDW_RADII_MEDIAN
        warnings.warn(f"Van-der-Waals radius for {element} unknown. Using median {radius:.2f}.")
    return radius
