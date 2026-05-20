# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Custom error types."""


class LowCoordinationNumber(KeyError):
    """Error for low coordination number."""


class HighCoordinationNumber(KeyError):
    """Error for high coordination number."""


class NoOpenDefined(KeyError):
    """Error in case the open check is not defined for this coordination number."""
