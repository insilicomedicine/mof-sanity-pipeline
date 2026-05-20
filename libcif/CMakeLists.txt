# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
cmake_minimum_required(VERSION 3.22)

project(libcif VERSION 1.0.0
               DESCRIPTION "Operations on cif files"
               LANGUAGES CXX)

# Temporary solution without presets
# Set optimization flag
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------

add_library(libcif STATIC)

target_sources(libcif PRIVATE src/cif.cpp
                              src/cif_ops.cpp
                              src/graph_ops.cpp
                              src/transform_ops.cpp
                              src/build_ops.cpp)

target_include_directories(libcif PUBLIC include/libcif)

target_compile_features(libcif PUBLIC cxx_std_20)

# ---------------------------------------
# ---- External dependencies options ----
# ---------------------------------------

option(LIBCIF_OBABEL_SMILES "Use OpenBabel for xyz to smiles conversion")
option(LIBCIF_RDK_CANON "Use RDKit for canonical ordering calculations")

# ---- openbabel support ----

if(LIBCIF_OBABEL)
  find_package(OpenBabel3 REQUIRED)

  add_compile_definitions(LIBCIF_BUILD_OBABEL)
  target_sources(libcif PRIVATE src/obabel_wrappers.cpp)
  target_link_libraries(libcif ${OpenBabel3_LIBRARIES})
  target_include_directories(libcif PUBLIC ${OpenBabel3_INCLUDE_DIRS})
endif()

# ---- rdkit support ----

if(LIBCIF_RDK)
  find_package(RDKit REQUIRED)
  find_package(Boost COMPONENTS timer system serialization iostreams REQUIRED)

  add_compile_definitions(LIBCIF_BUILD_RDKIT)
  target_sources(libcif PRIVATE src/mutation_ops.cpp src/rdkit_wrappers.cpp)
  target_link_libraries(libcif RDKit::SmilesParse RDKit::FileParsers RDKit::DetermineBonds)
  target_include_directories(libcif PUBLIC ${RDKit_INCLUDE_DIRS})
endif()

if(LIBCIF_SPGLIB)
  find_package(Spglib REQUIRED)

  add_compile_definitions(LIBCIF_BUILD_SPGLIB)
  target_sources(libcif PRIVATE src/spglib_wrappers.cpp)
  target_link_libraries(libcif symspg)
  target_include_directories(libcif PUBLIC ${Spglib_INCLUDE_DIRS})
endif()

if(LIBCIF_BLAKE2)
  add_compile_definitions(LIBCIF_BUILD_HASH)
  add_compile_definitions(USE_SSE)
  set_source_files_properties(external/sse/blake2b.c PROPERTIES LANGUAGE CXX)
  target_sources(libcif PRIVATE external/sse/blake2b.c external/compute_blake2b_hash.cpp)
  target_include_directories(libcif PRIVATE external)
endif()

# ------------------
# ---- Examples ----
# ------------------

add_subdirectory(examples/)
