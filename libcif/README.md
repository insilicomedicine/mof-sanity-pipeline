# libcif: A C++ Library for MOF Analysis and Design

## 1. Overview

**libcif** is a C++ library for the analysis, modification, and construction of Metal-Organic Frameworks (MOFs). It provides a toolkit for computational materials science applications. The library parses Crystallographic Information Files (CIFs), represents their data in internal structures, and offers functionality ranging from property calculation to structural modification.

The library is built with C++20 and can interface with external toolkits like RDKit, OpenBabel, and Spglib to extend its capabilities.

## 2. Core Functionality

### 2.1. CIF Parsing and Data Handling

The library includes functionality for parsing CIF files. It reads a `.cif` file and separates its contents into `data` and `loop` blocks. This information is used to populate a `cif_data` object, which stores:

*   **Unit Cell Parameters:** Lattice constants (`a`, `b`, `c`) and angles (`alpha`, `beta`, `gamma`).
*   **Atomic Coordinates:** The library handles both fractional and Cartesian (`Cartn`) coordinates for all atoms in the unit cell.
*   **Supercell Generation:** It can generate a supercell of arbitrary size (e.g., 3x3x3) to analyze bonds that cross unit cell boundaries.

### 2.2. Structural Analysis and Validation

`libcif` can build a molecular graph representation of a MOF to analyze its connectivity and chemical validity.

*   **Graph Construction:** A graph is built where atoms are nodes and chemical bonds are edges. The library can use different strategies (e.g., based on Cordero covalent radii) and a user-defined tolerance to determine connectivity.
*   **Structural Integrity Checks:** The `cif_is_good_v2` function validates a structure by checking for:
    *   **Overlapping Atoms:** Ensures no two atoms are closer than a minimum distance (`r_min`).
    *   **Disconnected Components:** Verifies that the framework is a single, connected component.
    *   **Correct Valency:** Checks that atoms do not exceed their maximum allowed number of bonds.
    *   **Geometric Soundness:** Performs checks on bond angles to validate the hybridization (sp, sp2, sp3) of carbon and nitrogen atoms.

### 2.3. Property Calculation

The library provides functions to calculate fundamental physical properties of a MOF directly from its CIF data:

*   **Chemical Formula:** Calculates the stoichiometric formula of the unit cell.
*   **Cell Volume:** Computes the volume of the unit cell in cubic angstroms.
*   **Mass and Density:** Calculates the molar mass of the unit cell and the crystal's density in g/cm³.

### 2.4. MOF Decomposition into Building Blocks

`libcif` includes functionality to decompose a MOF structure into its building blocks: inorganic **nodes** and organic **linkers**.

The `decompose_mof_v3` function performs this by:
1.  Expanding the unit cell into a supercell.
2.  Building a complete molecular graph of the supercell.
3.  Identifying all metal atoms.
4.  Using a Depth-First Search (DFS) algorithm to traverse the graph, separating atoms belonging to linkers from those belonging to nodes.
5.  Identifying and separating any disconnected fragments (e.g., solvent molecules).
6.  Using a built-in Blake2b hashing algorithm (`LIBCIF_BLAKE2`) to find and count the number of unique nodes and linkers.
7.  The final output consists of `.xyz` files for each unique building block.

### 2.5. *De Novo* MOF Construction

The library also provides tools to construct new MOF structures using a building-block approach. The functions in `build_ops.cpp` allow a user to:
1.  Define nodes and linkers as blueprint `.xyz` files.
2.  Specify a target topology that defines how these blocks should be connected.
3.  The library then uses geometric transformations (rotations calculated from Euler angles and translations) to assemble the building blocks into the final crystal structure.

### 2.6. Chemical Mutation Framework

`libcif` includes a framework for performing targeted chemical modifications on MOFs. This system relies on the RDKit library (`LIBCIF_RDK`) and uses a "gene" concept to represent potential mutation sites. It can:

*   **Identify Equivalent Sites:** By using RDKit's canonical atom ranking, it can identify chemically equivalent atoms within a MOF.
*   **Target Functional Groups:** The library can find specific functional groups within aromatic rings, such as C-H, C-F, or -CH3 groups.
*   **Execute Mutations:** It allows for systematic chemical transformations, such as substituting a metal ion or replacing an atom in a ring.

## 3. Building the Library

`libcif` uses CMake for its build system. The build process allows you to enable features based on which external libraries are available on your system.

### 3.1. Prerequisites

*   A C++20 compatible compiler (e.g., GCC 10+, Clang 12+).
*   CMake (version 3.22 or higher).
*   Optional external libraries (see below).

### 3.2. Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repository_url>
    cd libcif
    ```

2.  **Configure with CMake:**
    Create a build directory and run CMake. Use `-D` flags to enable optional features.
    ```bash
    # Basic build
    cmake -S . -B build --fresh

    # Build with all features enabled
    cmake -S . -B build --fresh -DLIBCIF_OBABEL=ON -DLIBCIF_RDK=ON -DLIBCIF_SPGLIB=ON -DLIBCIF_BLAKE2=ON
    ```

3.  **Compile the code:**
    ```bash
    cmake --build build -j
    ```
    The library (`libcif.a`) and example executables will be located in the `build/` directory.

### 3.3. CMake Build Flags

The following flags can be passed to CMake to enable optional functionality:

| Flag                 | Required Library | Description                                                                                                                                      |
| -------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `LIBCIF_OBABEL`      | OpenBabel 3      | Enables conversion of atomic structures to canonical SMILES strings. Required for the `cif_to_linker` example.                                   |
| `LIBCIF_RDK`         | RDKit, Boost     | Enables the chemical mutation framework and canonical atom ordering.                                                                             |
| `LIBCIF_SPGLIB`      | Spglib           | Enables crystallographic operations, such as finding the primitive unit cell of a structure.                                                      |
| `LIBCIF_BLAKE2`      | (none)           | Compiles a built-in Blake2b hashing implementation. Required for the MOF decomposition workflow (`cif_to_building_blocks_v3`) to identify unique components. |

## 4. Usage Examples

See [examples/README.md](./examples/README.md)

## 5. License and Citing

`libcif` is distributed as part of the MOF Sanity Pipeline and is released
under the same MIT License as the rest of the repository. License terms,
third-party attributions (BLAKE2, OpenBabel, RDKit, Spglib, Boost), and
citation instructions are documented in the root [`README.md`](../README.md)
and [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md).

