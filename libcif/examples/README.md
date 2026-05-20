# Examples description

The following examples demonstrate how to use
  the compiled executables to perform common tasks.
The examples are divided into 4 categories of their functionality.

Some more interesting examples, including linker mutations and others,
  can be found in `examples/legacy/` folder.
They are not documented here.

## Common parameters for graph calculation

The following parameters are common between many executables, because
  they are related to molecular graph calculation.
Here is the brief overview, later in the text they will be
  referred to as following.

* `r_min`: The minimum allowed distance between any two atoms.
* `r_max`: The maximum allowed bond length.
* `tolerance`: A tolerance factor for determining covalent bonds. See below for explanation.
* `strategy`: The covalent radii dictionary to use.

Explanations:

1. Notice, that `r_min` and `r_max` parameters are not used in all executables
  they're present in,
  and at the moment are kept for legacy and backwards compatibility reasons.
2. Tolerance parameters works as follows. If given a *positive* real number,
  two atoms are considered to be bonded, **IF** `r < tolerance*(cov_rad_1 + cov_rad_2)`.
If given a *negative* number, then bonded
  **IF** `r < |tolerance| + (cov_rad_1 + cov_rad_2)`
3. Strategy might be one of the {*cordero*, *pykko*, *jmol*}.

**RECOMMENDED DEFAULTS:**

```
r_min r_max tolerance strategy = 0.5 5.0 -0.45 jmol
```

or

```
r_min r_max tolerance strategy = 0.5 5.0 1.3 (defaults to cordero)
```

---

## 1. Validate a CIF Structure

The following tools check MOFs for structural and chemical validity.

### cif_ok

* **Usage:**

    ```shell
    ./cif_ok_v2 <filename.cif> <r_min> <r_max> <tolerance> [strategy]
    ```

* **Output:**
  * *stdout*: OK/BAD
  * *stderr:* Error messages

### cif_ok_v2

Same as `cif_ok`, but additional experimental checks for
  carbon and nitrogen geometries

### cif_check

* **Usage:** same as `cif_ok`
* **Output:**
  * *stdout:*
  `OK/BAD brute_formula cell_volume(angstrom^3)
    density(g/cm^3) a b c alpha beta gamma`
  * *stderr:* Error messages

### cif_check_v2

Same as `cif_check`, but uses `cif_ok_v2` under the hood.

### cif_quality

* **Usage:**

    ```shell
    ./cif_ok_v2 <filename.cif> <r_min> <r_max> <tolerance> <strategy> <penalty>
    ```

* **Output:**
  * *stdout*: (n_bad_atoms/n_atoms)^penalty
  * *stderr:* Error messages


## 2. Decompose a MOF into Building Blocks

These tools decompose a given MOF into its constituent linkers and nodes.

### cif_to_building_blocks_v2

* **Usage:**

    ```bash
    ./cif_to_building_blocks_v3 <filename.cif> <r_min> <r_max> <tolerance> [strategy]
    ```

* **Output:** This tool generates several output files in the working directory:
  * `linker_N.xyz`: An XYZ file for each unique linker found.
  * `node_N.xyz`: An XYZ file for each unique node found.
  * `disconnected_N.xyz`: An XYZ file for each disconnected component found.
  * `num_building_blocks.txt`: A summary of the number of each
  * `linker_N_infinite.xyz` (node, disconnected + num_building_blocks):
    An XYZ files for infinite components (if any).
    unique building block in the structure.
  * *stderr:* error messages

### cif_to_building_blocks_v3

Same as `cif_to_building_blocks_v2`, but uses linear-complexity graph calculation

## 3. Calculate different properties utilities

### calc_charge

Calculates the total charge of a cif file (if charges field is given, otherwise zero).

* **Usage:** `./calc_charge <filename.cif>`
* **Output:**
  * *stdout:* total charge of the cif
  * *stderr:* errors

### cif_hash

Calculates the hash of a cif structure in different ways.
**Requires `-DLIBCIF_BLAKE2=ON` during configuration.**

* **Usage:** `./cif_hash <filename.cif>`
* **Output:**
  * *stdout:* hashes of the cif
  * *stderr:* errors

### cif_prim_spg

Calculates the primitivized unit cell.
**Requires `-DLIBCIF_SPGLIB=ON` during configuration.

* **Usage:** `./cif_prim_spg <filename.cif>`
* **Output:**
  * *stdout:* new cif file
  * *stderr:* errors

## 4. Build new MOFs with specific geometries from SBUs

The following utilities build new MOFs given the new xyz file for linkers.
They are specifically tailored to build particular structures, so they might not
  be generalizable, but illustrate additions to the library.
They utilize the transformations module extensively, specifically
  the Kabsch algorithm to align linkers to the destination.

1. build_but.cpp
2. build_fex.cpp
3. make_999.cpp
4. make_sifsix.cpp
5. build_material-10.cpp
