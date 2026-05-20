# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.

"""Version information for :mod:`framechecker`.

Run with ``python -m framechecker.version``
"""

import os
from subprocess import CalledProcessError, check_output  # noqa: S404

__all__ = [
    "VERSION",
    "get_version",
    "get_git_hash",
]

VERSION = "0.9.7-dev"


def get_git_hash() -> str:
    """Get the :mod:`sonipredict` git hash."""
    with open(os.devnull, "w") as devnull:
        try:
            ret = check_output(  # noqa: S603,S607
                ["git", "rev-parse", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=devnull,
            )
        except CalledProcessError:
            return "UNHASHED"
        else:
            return ret.strip().decode("utf-8")[:8]


def get_version(with_git_hash: bool = False):
    """Get the :mod:`sonipredict` version string, including a git hash."""
    return f"{VERSION}-{get_git_hash()}" if with_git_hash else VERSION


if __name__ == "__main__":
    print(get_version(with_git_hash=True))  # noqa:T201
