# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Custom error types."""


class NoMetal(KeyError):
    """Error in case there is no metal in structure."""
