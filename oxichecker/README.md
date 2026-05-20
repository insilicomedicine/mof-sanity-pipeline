# OxiChecker

OxiChecker is the Stage-2 charge-balance and decomposition validator of the MOF
Sanity Pipeline. It takes a CIF, decomposes it into nodes and linkers via
`libcif`'s `cif_to_building_blocks_v3` binary, assigns formal oxidation
states / fragment charges, and reports whether the structure is charge-balanced
(and whether any of its building blocks are flagged as infinite).

For the role OxiChecker plays inside the full pipeline, the per-row meanings of
the `OxiChecker Validity` column, and how its verdict combines with the other
stages into the final `Sanity` value, see the root [`README.md`](../README.md).

## Layout

```
oxichecker/
├── main.py              # CLI entrypoint; orchestrates per-structure validation
├── fragment_analyzer.py # node fragment charge assignment (water, hydroxide,
│                        # nitrates, halides, oxos, peroxides, ...)
├── linker_analyzer.py   # linker charge assignment from .xyz (COO, PO3, NH2,
│                        # pillar metals, perchlorates, ...)
└── entrypoint.sh        # in-container UID/GID resolver used when OxiChecker
                         # is invoked as a standalone container
```

`main.py` calls `cif_to_building_blocks_v3` from `libcif` to obtain per-CIF
`.xyz` files for nodes and linkers, then feeds those `.xyz` files through
`fragment_analyzer` and `linker_analyzer` to compute the total charge of the
unit cell. A structure passes when the summed charge equals zero (with the
"infinite node / linker" caveats listed in the root README).

## Standalone usage

OxiChecker is normally driven by `sanity_runner` and the full pipeline CLI
described in the root README. It can also be invoked directly on a folder of
P1-symmetry CIFs:

```bash
python main.py --input <CIFS_DIR> --output result.csv --num_cores <N>
```

### Key flags

| Flag | Default | Meaning |
|---|---|---|
| `--input` | *(required)* | CIF file or directory of CIFs |
| `--output` | `result.csv` | Output CSV path |
| `--version` | `v3` | Decomposer to use: `v2`, `v3`, or `external` |
| `--cif-to-bb-path` | — | External `cif_to_building_blocks` binary (used with `--version external`) |
| `--num_cores` | auto | Parallel worker count |
| `--nosave-decompose-data` | off | Delete per-structure decomposition artefacts after the run |
| `--decompose-dir` | `./decompose_<output_name>` | Where to keep intermediate `.xyz` files |
| `--verbose` | off | Log every decomposition command executed |
| `--param1` ... `--param4` | `0.5`, `5.0`, `-0.45`, `jmol` | Decomposer tuning parameters |
| `--general-timeout` | `300` | Per-structure cap (seconds) |
| `--decomposition-timeout` | `60` | Cap for the `cif_to_building_blocks_v3` call |
| `--obabel-timeout` | `30` | Cap for the OpenBabel conversion step |

When OxiChecker is invoked through `sanity_runner`, these flags are exposed
with an `--oxichecker-` prefix (e.g. `--oxichecker-num_cores`). See the root
README for the full pipeline-level flag list.

## License and Citing

OxiChecker is part of the MOF Sanity Pipeline and is distributed under the
same MIT License as the rest of the repository. License terms, third-party
attributions, and citation instructions are documented in the root
[`README.md`](../README.md) and [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md).
