# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Types reused across the package."""
from pathlib import Path
from typing import Union

from pymatgen.core import IStructure, Structure
from typing_extensions import TypeAlias

PathType: TypeAlias = Union[str, Path]
StructureIStructureType: TypeAlias = Union[Structure, IStructure]
