# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# File from MOFChecker (https://github.com/lamalab-org/mofchecker).
# Copyright (c) 2021 Kevin Jablonka -- MIT License.
# See repository LICENSE and THIRD_PARTY_NOTICES.md for full terms.
"""Command line interface."""

import json

import click

from framechecker import DESCRIPTORS, FrameChecker


@click.command()
@click.option(
    "--primitive/--no-primitive",
    default=True,
    help="Perform the analysis on the primitive structure",
    show_default=True,
)
@click.option(
    "--descriptors",
    "-d",
    multiple=True,
    type=click.Choice(DESCRIPTORS),
    default=DESCRIPTORS,
    help="Select descriptors to be computed.",
    show_default=False,
)
@click.argument("CIF_FILES", type=click.Path(exists=True, dir_okay=False), nargs=-1)
def run(primitive, descriptors, cif_files):
    """Check provided structures and print list of JSON objects with descriptors."""
    print("[")  # noqa: T201
    for index, structure_file in enumerate(cif_files):
        framechecker = FrameChecker.from_cif(structure_file, primitive=primitive)
        descriptors = framechecker.get_frame_descriptors(descriptors=descriptors)

        string = json.dumps(descriptors, indent=2)
        if index != len(cif_files) - 1:
            string += ","
        print(string)  # noqa: T201
    print("]")  # noqa: T201
