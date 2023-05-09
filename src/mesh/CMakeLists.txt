# Copyright (c) 2019, Triad National Security, LLC
# All rights reserved.

# Copyright 2019. Triad National Security, LLC. This software was
# produced under U.S. Government contract 89233218CNA000001 for Los
# Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S. Department of Energy. 
# All rights in the program are reserved by Triad National Security,
# LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting
# on its behalf a nonexclusive, paid-up, irrevocable worldwide license
# in this material to reproduce, prepare derivative works, distribute
# copies to the public, perform publicly and display publicly, and to
# permit others to do so
 
# 
# This is open source software distributed under the 3-clause BSD license.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of Triad National Security, LLC, Los Alamos
#    National Laboratory, LANL, the U.S. Government, nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
#  
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#
#  Jali
#    Mesh 
#
#

# Jali module, include files found in JALI_MODULE_PATH
# include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(JALI_MESH)

set(JALI_MESH_headers
  MeshDefs.hh
  Mesh.hh
  MeshTile.hh
  MeshSet.hh
  block_partition.hh
  )
list(TRANSFORM JALI_MESH_headers PREPEND "${JALI_MESH_SOURCE_DIR}/")

set(JALI_MESH_sources
  Mesh.cc
  MeshTile.cc
  MeshSet.cc
  block_partition.cc
  )


add_library(jali_mesh ${JALI_MESH_sources})
set_target_properties(jali_mesh PROPERTIES PUBLIC_HEADER "${JALI_MESH_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_mesh ALIAS jali_mesh)


# We should change install_interface to include/jali/jali_mesh
# but that requires changing all client codes
target_include_directories(jali_mesh PUBLIC
  $<BUILD_INTERFACE:${JALI_MESH_BINARY_DIR}>
  $<BUILD_INTERFACE:${JALI_MESH_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )


#
# Library: mesh
#

target_link_libraries(jali_mesh PUBLIC jali_geometry)
target_link_libraries(jali_mesh PUBLIC jali_error_handling)


# Factory class
add_subdirectory(mesh_factory)

# Derived and support classes

add_subdirectory(mesh_simple)
target_link_libraries(jali_mesh PUBLIC jali_simple_mesh)

# Mesh Frameworks

# STK (Trilinos Package)
#if (ENABLE_STK_Mesh)
#    add_subdirectory(mesh_stk)
#
#  target_link_libraries(jali_mesh PUBLIC jali_stk_mesh)
#  target_compile_defintions(jali_mesh PUBLIC HAVE_STK_MESH)
#endif(ENABLE_STK_Mesh)    

# MOAB Mesh framework
#if (ENABLE_MOAB_Mesh)
#    add_subdirectory(mesh_moab)
#
#    target_link_libraries(jali_mesh PUBLIC jali_moab_mesh)
#    target_compile_defintions(jali_mesh PUBLIC HAVE_MOAB_MESH)
#endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
  add_subdirectory(mesh_mstk)

  target_link_libraries(jali_mesh PUBLIC jali_mstk_mesh)
  target_compile_definitions(jali_mesh PUBLIC HAVE_MSTK_MESH)
endif(ENABLE_MSTK_Mesh)    





# Add third party libraries as dependencies (found in JaliTPL.cmake
# through top level CMakeList.txt)
if (NOT METIS_LIBRARIES OR NOT TARGET ${METIS_LIBRARIES})
  find_package(METIS QUIET REQUIRED)   # uses MSTK module
  message(STATUS "Found Metis library: ${METIS_LIBRARY}")
endif ()

# Make METIS a dependency of jali_mesh
target_compile_definitions(jali_mesh PUBLIC Jali_HAVE_METIS)
target_link_libraries(jali_mesh PUBLIC ${METIS_LIBRARIES})

if (ENABLE_ZOLTAN)
  if (NOT Zoltan_LIBRARIES OR NOT TARGET ${Zoltan_LIBRARIES})
    # Starting from CMake 3.15 we can tell find_package we prefer for it
    # to look for config files first and then fall back to modules
    find_package(Zoltan QUIET REQUIRED)
    message(STATUS "Found TARGET Zoltan_LIBRARIES = ${Zoltan_LIBRARIES}")
  endif()

  # Make Zoltan a dependency of jali_mesh
  target_compile_definitions(jali_mesh PUBLIC Jali_HAVE_ZOLTAN)
  target_link_libraries(jali_mesh PUBLIC ${Zoltan_LIBRARIES})

  # Zoltan's cmake config does not write out include dir dependency for
  # the Zoltan library. So explicitly enable it
  target_include_directories(jali_mesh PUBLIC ${Zoltan_INCLUDE_DIRS})
endif (ENABLE_ZOLTAN)

if (NOT ExodusII_LIBRARIES OR NOT TARGET ${ExodusII_LIBRARIES})
  # First seee if a config file got installed as part of the SEACAS project
  # NOTE: NOT ABLE TO PROCESS THIS CORRECTLY
  #  find_package(SEACASExodus QUIET CONFIG PATHS ${ExodusII_DIR})
  
  #  if (SEACASExodus_FOUND)
  #    target_include_directories(${MSTKLIB} PUBLIC ${ExodusII_INCLUDE_DIRS})
  #    target_link_libraries(${MSTKLIB} PUBLIC ${ExodusII_LIBRARIES})
  #    message(STATUS "Found Exodus II library (from SEACAS): ${SEACAS_LIBRARIES}")
  #  else ()
  find_package(ExodusII QUIET REQUIRED)  # Uses custom module
  message(STATUS "Found Exodus II library: ${ExodusII_LIBRARY}")
endif ()

# Make ExodusII a dependency of jali_mesh (strictly speaking it is not...)
target_link_libraries(jali_mesh PUBLIC ${ExodusII_LIBRARIES})


install(TARGETS jali_mesh
  EXPORT JaliTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

if (BUILD_TESTS)
  
  # Need to copy files for the tests 
  if (NOT (${JALI_MESH_SOURCE_DIR} EQUAL ${JALI_MESH_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${JALI_MESH_SOURCE_DIR}/test ${JALI_MESH_BINARY_DIR}/test) 
  endif()

  # Test: mesh geometry

  add_Jali_test(mesh_geometry test_mesh_geometry
    KIND unit
    SOURCE test/Main.cc test/test_mesh_geometry.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  # Test mesh entities

  add_Jali_test(mesh_iterators test_entity_iterators
    KIND unit
    SOURCE test/Main.cc test/test_entity_iterators.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  add_Jali_test(mesh_entity_counts_serial mesh_entity_counts_serial
    KIND unit
    SOURCE test/Main.cc test/test_entity_counts.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  add_Jali_test(mesh_entity_counts_parallel mesh_entity_counts_parallel
    KIND unit
    NPROCS 4
    SOURCE test/Main.cc test/test_entity_counts.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  # Test side entities
  
  add_Jali_test(side_tests test_sides
    KIND unit
    SOURCE test/Main.cc test/test_sides.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  # Test wedge entities
  
  add_Jali_test(wedge_tests test_wedges
    KIND unit
    SOURCE test/Main.cc test/test_wedges.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  # Test corner entities
  
  add_Jali_test(corner_tests test_corners
    KIND unit
    SOURCE test/Main.cc test/test_corners.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  # Test mesh tiles
  
  add_Jali_test(tile_tests test_one_tile
    KIND unit
    SOURCE test/Main.cc test/test_one_tile.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  add_Jali_test(mesh_tile_tests_serial test_meshtiles_serial
    KIND unit
    SOURCE test/Main.cc test/test_meshtiles.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  add_Jali_test(mesh_tile_tests_parallel test_meshtiles_parallel
    KIND unit
    NPROCS 4
    SOURCE test/Main.cc test/test_meshtiles.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})
  
  # Test boundary ghosts
  
  add_Jali_test(mesh_boundary_ghost_tests_serial test_boundary_ghosts_serial
    KIND unit
    SOURCE test/Main.cc test/test_boundary_ghosts.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  add_Jali_test(mesh_boundary_ghost_tests_parallel test_boundary_ghosts_parallel
    KIND unit
    NPROCS 4
    SOURCE test/Main.cc test/test_boundary_ghosts.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})


  # Test mesh sets
  
  add_Jali_test(mesh_set_tests_serial test_meshsets_serial
    KIND unit
    SOURCE test/Main.cc test/test_meshsets.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  add_Jali_test(mesh_set_tests_parallel test_meshsets_parallel
    KIND unit
    NPROCS 4
    SOURCE test/Main.cc test/test_meshsets.cc
    LINK_LIBS jali_mesh jali_mesh_factory ${UnitTest++_LIBRARIES})

  # test block partitioning
  add_Jali_test(block_partitioning block_partitioning
    KIND unit
    SOURCE test/Main.cc test/test_block_partition.cc
    LINK_LIBS jali_mesh ${UnitTest++_LIBRARIES})

endif()
  
