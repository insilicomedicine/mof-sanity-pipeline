# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Helper functions for the check functions."""
from shutil import which


def is_tool(name: str) -> bool:
    """Check whether `name` is on PATH and marked as executable.

    Based onhttps://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script

    Args:
        name (str): The name of the tool to check.

    Returns:
        bool: Whether the tool is on PATH.
    """
    return which(name) is not None
