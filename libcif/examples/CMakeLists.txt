# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
cmake_minimum_required(VERSION 3.22)

project(libcif_examples)

# *************
# LEGACY examples are commented out
# *************

add_executable(cif_ok cif_ok.cpp)
add_executable(cif_ok_v2 cif_ok_v2.cpp)
add_executable(cif_quality cif_quality.cpp)
add_executable(cif_check cif_check.cpp)
add_executable(cif_check_v2 cif_check_v2.cpp)
add_executable(calc_charge calc_charge.cpp)
target_link_libraries(cif_ok PRIVATE libcif)
target_link_libraries(cif_ok_v2 PRIVATE libcif)
target_link_libraries(cif_check PRIVATE libcif)
target_link_libraries(cif_check_v2 PRIVATE libcif)
target_link_libraries(cif_quality PRIVATE libcif)
target_link_libraries(calc_charge PRIVATE libcif)

# Material conctruction
add_executable(build_uit build_uit.cpp)
add_executable(build_fex build_fex.cpp)
add_executable(build_but build_but.cpp)
add_executable(build_sifsix build_sifsix.cpp)
target_link_libraries(build_uit PRIVATE libcif)
target_link_libraries(build_fex PRIVATE libcif)
target_link_libraries(build_but PRIVATE libcif)
target_link_libraries(build_sifsix PRIVATE libcif)

# add_executable(cif_ok_multi_test cif_ok_multi_test.cpp)
# target_link_libraries(cif_ok_multi_test PRIVATE libcif)
# add_executable(tmp_material-10 tmp_material-10.cpp)
# target_link_libraries(tmp_material-10 PRIVATE libcif)
# add_executable(cif_to_linker_xyz cif_to_linker_xyz.cpp)
# target_link_libraries(cif_to_linker_xyz PRIVATE libcif)
#

if(LIBCIF_SPGLIB)
  add_executable(cif_prim_spg cif_prim_spg.cpp)
  target_link_libraries(cif_prim_spg PRIVATE libcif)
endif()

if(LIBCIF_BLAKE2)
  add_executable(cif_hash cif_hash.cpp)
  add_executable(cif_to_building_blocks_v2 cif_to_building_blocks_v2.cpp)
  add_executable(cif_to_building_blocks_v3 cif_to_building_blocks_v3.cpp)
  target_link_libraries(cif_hash PRIVATE libcif)
  target_link_libraries(cif_to_building_blocks_v2 PRIVATE libcif)
  target_link_libraries(cif_to_building_blocks_v3 PRIVATE libcif)

  # add_executable(cif_hash_multi cif_hash_multi.cpp)
  # target_link_libraries(cif_hash_multi PRIVATE libcif)
endif()

# LEGACY only down
# Still might be interesting, especially the RDKIT-assisted mutation
# part. Some work has been done.

if(LIBCIF_OBABEL)
  add_executable(cif_to_linker cif_to_linker.cpp)
  add_executable(cif_to_linker_v2 cif_to_linker_v2.cpp)
  target_link_libraries(cif_to_linker PRIVATE libcif)
  target_link_libraries(cif_to_linker_v2 PRIVATE libcif)
endif()

if(LIBCIF_RDK)
  add_executable(make_chromosome make_chromosome.cpp)
  add_executable(canon_order canon_order.cpp)
  add_executable(mutate mutate.cpp)
  add_executable(determine_charge determine_charge.cpp)
  target_link_libraries(make_chromosome PRIVATE libcif)
  target_link_libraries(canon_order PRIVATE libcif)
  target_link_libraries(mutate PRIVATE libcif)
  target_link_libraries(determine_charge PRIVATE libcif)
endif()

