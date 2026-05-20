# MOF Sanity Pipeline

Docker-based CIF sanity checker pipeline combining **LibCif**, **PLATON**, **MOFChecker** and **OxiChecker** for validation, decomposition and charge-balance analysis of metal-organic frameworks (MOFs).

The container entrypoint provides a unified CLI that orchestrates all tools from a single command.

---

## Requirements

- Docker
- PLATON source files (see below — they must be downloaded separately due to licensing)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/insilicomedicine/mof-sanity-pipeline.git
cd mof-sanity-pipeline
```

### 2. Obtaining PLATON sources

PLATON is **not distributed** with this repository due to its license terms.
Obtaining the sources is the user's responsibility: visit the official PLATON
site, review and agree to the applicable license, and then download the
source files. PLATON is free of charge for academic use; for-profit
organisations must contact [Ton Spek](mailto:a.l.spek@uu.nl) (or submit the
PLATON application form referenced on the website) and obtain a license
before copying the sources.

Source: **https://www.platonsoft.nl/platon/pl030000.html**

Once licensing is settled, download the two files below and place them in the
`code/` directory:

- `platon.f.gz`
- `xdrvr.c.gz`

### 3. Build the Docker image

```bash
docker build . -t sanity-pipeline
```

The container resolves UID/GID at runtime by inspecting the mounted output directory, so no host UID needs to be baked into the image. Build once, run as any user.

Optional build-time flags:

- `--build-arg RUN_ENTRYPOINT_TESTS=yes` — run CLI smoke tests during build (default: no)

---

## Usage

Mount your working directory into `/data` and point the entrypoint at the input and output subpaths inside it. The container resolves UID/GID from the owner of `/data` and drops privileges to that user via `gosu` before any work is done, so the output directory is created from inside the container with your UID/GID — no host-side `mkdir` and no root-owned artefacts:

```bash
docker run --rm \
    -v "$PWD":/data \
    sanity-pipeline \
    -i /data/<CIF_FOLDER> \
    --output-dir /data/my_results \
    --n-jobs <NUMBER_OF_CPUS>
```

`./my_results/` appears in `$PWD` after the run, owned by you. If `--output-dir` is omitted, results land in `<input>_sanity_results/` next to the input (same ownership rules apply).

### Run modes

| Flag | Meaning |
|---|---|
| *(no flag)* | Full pipeline: sanity runner → OxiChecker → postprocessing |
| `--run-sanity-only` | LibCif + PLATON + MOFChecker validation only |
| `--run-oxichecker-only` | Charge balance / oxidation state validation only |
| `--run-postprocess-only` | Final CSV consolidation only |

### Common flags

| Flag | Meaning |
|---|---|
| `-i`, `--input VALUE` | CIF file or directory |
| `--n-jobs VALUE` | Parallel worker count (forwarded to all stages) |
| `--output-dir VALUE` | Custom output directory |

### Forwarding flags to underlying tools

Use prefixes to forward options to each stage:

- `--sanity-*` for sanity runner (e.g. `--sanity-total-timeout 600`)
- `--oxichecker-*` for OxiChecker (e.g. `--oxichecker-version v3`)
- `--post-*` for postprocessing (e.g. `--post-output my.csv`)

See `docker run --rm sanity-pipeline --help` for the full flag list.

---

## Output Structure

Main output directory: `{input_name}_sanity_results/`

```
{input_name}_sanity_results/
├── babel_cifs/                      # OpenBabel-standardised CIFs
├── decompose_oxichecker_results/    # OxiChecker decomposition per structure
├── p1_cifs/                         # P1-symmetry CIFs
├── results/                         # Per-structure validation JSONs and logs
├── symm_cifs/                       # Symmetry-processed CIFs
├── oxichecker_results.csv           # OxiChecker stage results
└── sanity_results.csv               # Final consolidated report
```

### `sanity_results.csv` — main columns

**Identification:** `cif`, `content_hash`, `formula`, `reduced_formula`

**Structural properties:** `density`, `volume`, `group_str`, `group_id`, `structure_hash_strict`, `structure_hash`

**Graph analysis:** `is_graph_constructed`, `graph_dim`, `HAS_OMS`, `DECORATED_GRAPH_HASH`, `UNDECORATED_GRAPH_HASH`, `DECORATED_SCAFFOLD_HASH`, `UNDECORATED_SCAFFOLD_HASH`

**Validation results:**
- `Basic Validity` — basic structural validation (boolean)
- `LibCif Validity` — LibCif fast filter (boolean)
- `LibCif_Warning` — LibCif warning messages
- `PLATON Validity` — PLATON geometry check (boolean)
- `OxiChecker Validity` — charge balance / oxidation state result (string, see below)
- `Sanity` — overall verdict (boolean)

**Chemical sanity checks** (`True` means problem detected):
`has atomic overlaps`, `has overcoordinated c/h/n`, `has lone molecule`, `has bad rare earth`, `has bad alkali alkaline`, `has bad terminal oxo`

### `OxiChecker Validity` values

| Value | Meaning |
|---|---|
| `True` | Successfully decomposed; charges balanced |
| `True; Infinity Node` | Valid; structure contains an infinite node |
| `True; Infinity Linker` | Valid; structure contains an infinite linker |
| `True; no metal` | Structure has no metals (e.g. COF); charge check skipped |
| `Invalid charges` | Charge balance failed |
| `Invalid charges; Infinity Node` | Charge balance failed; infinite node detected |
| `Invalid charges; Infinity Linker` | Charge balance failed; infinite linker detected |
| `Failed to decompose` | LibCif decomposer did not produce output |
| `Empty decomposition` | Decomposer returned no building blocks |
| `No building blocks file` | `num_building_blocks.txt` not found after decomposition |
| `Processing timeout` | Exceeded `--oxichecker-general-timeout` |
| `Processing error: ...` | Unhandled exception during processing |
| `Not processed` | OxiChecker stage was not run (e.g. with `--run-sanity-only`) |

### `LibCif Validity` values

| Value | Meaning |
|---|---|
| `True` | Passed LibCif filter |
| `False` | Failed LibCif filter |

### Overall `Sanity` verdict

A structure passes overall sanity when **all** conditions hold:

1. All chemical sanity checks are `False` (no problems detected)
2. `Basic Validity`, `LibCif Validity`, and `PLATON Validity` are not `False`
3. `OxiChecker Validity` either starts with `True` or equals `Not processed`

---

## Repository Structure

```
mof-sanity-pipeline/
|-- code/                       # Core Python pipeline (sanity_runner, postprocessing, utilities)
|-- oxichecker/                 # OxiChecker module (main.py, fragment_analyzer.py, linker_analyzer.py)
|-- libcif/                     # C++ LibCif library (compiled during build)
|-- scripts/                    # CLI smoke tests
|-- test/                       # Test CIF files
|-- Dockerfile
|-- .dockerignore               # excludes VCS, docs, caches and build artefacts from the build context
|-- docker-entrypoint.sh        # runtime UID/GID resolver
|-- generate_entrypoint.sh      # builds the in-container CLI dispatcher
|-- LLM_generated_cifs.zip      # 24,950 LLM-generated CIFs used in the paper (see Citing)
|-- THIRD_PARTY_NOTICES.md
|-- README.md
```

---

## LLM-generated CIFs

`LLM_generated_cifs.zip` contains 24,950 candidate MOF structures generated by a large language model. This is the exact dataset evaluated by the sanity pipeline in the paper referenced under [Citing](#citing); it is bundled here so the published results can be reproduced end-to-end. Unzip and pass the resulting `LLM_generated_cifs/` directory as the input to the pipeline:

```bash
unzip LLM_generated_cifs.zip
docker run --rm \
    -v "$PWD":/data \
    sanity-pipeline \
    -i /data/LLM_generated_cifs \
    --output-dir /data/llm_results \
    --n-jobs <NUMBER_OF_CPUS>
```

---

## Citing

If you use this pipeline, please cite the associated ChemRxiv preprint:

> Bezrukov, D., Pupeza, A., Younes, M., et al. *Sanity and Decomposition Pipeline for Metal-Organic Frameworks in Generative AI.* ChemRxiv (2026). DOI: [10.26434/chemrxiv.15003614/v1](https://doi.org/10.26434/chemrxiv.15003614/v1)

<details>
<summary>BibTeX</summary>

```bibtex
@article{Bezrukov2026,
  title = {Sanity and Decomposition Pipeline for Metal-Organic Frameworks in Generative AI},
  url = {http://dx.doi.org/10.26434/chemrxiv.15003614/v1},
  DOI = {10.26434/chemrxiv.15003614/v1},
  publisher = {American Chemical Society (ACS)},
  author = {Bezrukov, Dmitry and Pupeza, Aleksandr and Younes, Mourad and Rublev, Pavel and Romashin, Ivan and Kamorzin, Boris and Alhammad, Bashaer and Aljama, Hassan and Toro, Frankly and Alahmed, Ammar and Iamashev, Mikhail and Mazaleva, Olga and Permyakova, Anastasia and Aliper, Alex and Zhavoronkov, Alex and Badra, Jihad},
  year = {2026},
  month = May
}
```

</details>

## License

The code in this repository is released under the MIT License.

Parts of this code are derived from MOFChecker, which is also distributed
under the MIT License. The original MOFChecker copyright and license
notice are preserved in THIRD_PARTY_NOTICES.md.

PLATON is not included in this repository and is not redistributed under
this repository's MIT License. PLATON support is implemented only as an
integration layer that invokes the external PLATON binary. Users must
obtain PLATON separately and comply with the PLATON license terms.
