# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# Originally derived from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- original MOFChecker code, MIT License.
# Copyright (c) 2026 Saudi Aramco -- modifications, MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Tooling for finding open metal sites."""
from typing import List

import numpy as np
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import LocalStructOrderParams
from structuregraph_helpers.analysis import get_cn

from .definitions import OP_DEF
from .errors import HighCoordinationNumber, LowCoordinationNumber
from ..check_base import AbstractIndexCheck
from ..utils.get_indices import get_metal_indices
from ...errors import NoMetal
from ...types import StructureIStructureType


class FrameOMS(AbstractIndexCheck):
    """A 'checker' for finding open metal sites."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the FrameOMS class.

        Args:
            structure (StructureIStructureType): pymatgen structure
            structure_graph (StructureGraph): pymatgen structure graph
        """
        self.structure = structure
        self.structure_graph = structure_graph
        self._metal_indices = get_metal_indices(structure)
        self._open_indices: set = set()
        self._has_oms = None
        self.metal_features = {}

    @property
    def name(self) -> str:
        """Return the name of the check."""
        return "OMS"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check if there are any open metal sites in the structure."

    def get_cn(self, index):
        """Return the coordination number."""
        return get_cn(self.structure_graph, index)

    @classmethod
    def from_framechecker(cls, framechecker):
        """Initialize a OMS check from a framechecker instance."""
        checker = cls(framechecker.structure, framechecker.graph)
        checker.get_cn = framechecker.get_cn
        return checker

    def get_metal_descriptors_for_site(self, site_index: int) -> dict:
        """Compute the checks for one metal site."""
        if len(self._metal_indices) == 0:
            raise NoMetal
        return self._get_metal_descriptors_for_site(site_index)

    def _get_metal_descriptors(self):
        descriptordict = {}
        for site_index in self._metal_indices:
            descriptordict[site_index] = self._get_metal_descriptors_for_site(site_index)

        self.metal_features = descriptordict

        return descriptordict

    def get_metal_descriptors(self) -> dict:
        """Return local structure order parameters.

        Key is the site index.

        Raises:
            NoMetal: If no metal can be found in the structure

        Returns:
            dict: Key is the site index.
        """
        if len(self._metal_indices) == 0:
            raise NoMetal
        return self._get_metal_descriptors()

    def _run_check(self):
        indices = self.check_oms()
        return len(indices) == 0, indices

    def check_oms(self) -> List[int]:
        """Check if there are any open metal sites in the structure.

        True if the structure contains open metal sites (OMS).

        Also returns True in case of low coordination numbers (CN <=3)
        which typically indicate open coordination for Frames.
        For high coordination numbers, no good order parameter for open
        structures is available, and so we return `None` even though
        this might change in a future release.

        Raises:
            NoMetal: Raised if the structure contains no metal

        Returns:
            List[int]: OMS indices
        """
        oms_sites = []
#        if len(self._metal_indices) == 0:
#            raise NoMetal("This structure does not contain a metal")
        for site_index in self._metal_indices:
            if self.is_site_open(site_index):
                oms_sites.append(site_index)
        return oms_sites

    @staticmethod
    def _check_if_open(lsop, is_open, weights, threshold: float = 0.5):
        if lsop is not None:
            if is_open is None:
                return False
            lsop = np.array(lsop) * np.array(weights)
            open_contributions = lsop[is_open].sum()
            close_contributions = lsop.sum() - open_contributions
            return open_contributions / (open_contributions + close_contributions) > threshold
        return None

    def _get_metal_descriptors_for_site(self, site_index: int):
        metal = str(self.structure[site_index].species)
        try:
            (
                cn,  # pylint:disable=invalid-name
                names,
                lsop,
                is_open,
                weights,
            ) = self._get_ops_for_site(site_index)
            site_open = FrameOMS._check_if_open(lsop, is_open, weights)
            if site_open:
                self._open_indices.add(site_index)
            descriptors = {
                "metal": metal,
                "lsop": dict(zip(names, lsop)),
                "open": site_open,
                "cn": cn,
            }
        except LowCoordinationNumber:
            descriptors = {"metal": metal, "lsop": None, "open": True, "cn": None}
        except HighCoordinationNumber:
            descriptors = {"metal": metal, "lsop": None, "open": None, "cn": None}
        return descriptors

    def _get_ops_for_site(self, site_index):
        cn = self.get_cn(site_index)  # pylint:disable=invalid-name
        try:
            names = OP_DEF[cn]["names"]
            is_open = OP_DEF[cn]["open"]
            weights = OP_DEF[cn]["weights"]
            lsop = LocalStructOrderParams(names)
            return (
                cn,
                names,
                lsop.get_order_parameters(self.structure, site_index),
                is_open,
                weights,
            )
        except KeyError as exc:
            # For a bit more fine grained error messages
            if cn <= 3:  # pylint:disable=no-else-raise
                raise LowCoordinationNumber(
                    "Coordination number {} is low \
                        and order parameters undefined".format(
                        cn
                    )
                ) from exc
            elif cn > 8:
                raise HighCoordinationNumber(
                    "Coordination number {} is high \
                        and order parameters undefined".format(
                        cn
                    )
                ) from exc

            return cn, None, None, None, None

    def is_site_open(self, site_index: int) -> bool:
        """Check for a site if is open.

        This is based on the values of some coordination geometry fingerprints.

        Args:
            site_index (int): Index of the site in the structure

        Returns:
            bool: True if site is open
        """
        if site_index not in self._open_indices:
            try:
                _, _, lsop, is_open, weights = self._get_ops_for_site(site_index)
                site_open = FrameOMS._check_if_open(lsop, is_open, weights)
                if site_open:
                    self._open_indices.add(site_index)
                return site_open
            except LowCoordinationNumber:
                return True
            except HighCoordinationNumber:
                return None
        return True
