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
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH)

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})


# Check which partitioners are enabled and set some variables that
# can be checked in the code with an ifdef statement. ENABLE_METIS,
# ENABLE_ZOLTAN, METIS_DIR and ZOLTAN_DIR are all from Jali-tpl-config.cmake

if (ENABLE_METIS)
  message(STATUS "Adding definition Jali_HAVE_METIS")
  add_definitions(-DJali_HAVE_METIS)
  include_directories(${METIS_DIR}/include)
endif ()

if (ENABLE_ZOLTAN)
  message(STATUS "Adding definition Jali_HAVE_ZOLTAN")
  add_definitions(-DJali_HAVE_ZOLTAN)
  include_directories(${ZOLTAN_DIR}/include)
endif ()
if ((NOT ENABLE_METIS) AND (NOT ENABLE_ZOLTAN))
  message(STATUS "No partitioner enabled?")
endif()



#
# Library: mesh
#
set(mesh_link_libs geometry)
if (ENABLE_METIS)
  list(APPEND mesh_link_libs ${METIS_LIBRARIES})
endif (ENABLE_METIS)
if (ENABLE_ZOLTAN)
  list(APPEND mesh_link_libs ${ZOLTAN_LIBRARIES})
endif (ENABLE_ZOLTAN)

add_Jali_library(mesh SOURCE Mesh.cc MeshTile.cc MeshSet.cc 
  LINK_LIBS ${mesh_link_libs})

#
# Library: mesh_audit
#
#add_Jali_library(mesh_audit SOURCE MeshAudit.cc LINK_LIBS error_handling)

#
# Install Header files
#
file(GLOB mesh_inc_files "*.hh")
add_install_include_file(${mesh_inc_files})



# Derived and support classes

add_subdirectory(mesh_simple)

# Mesh Frameworks

# STK (Trilinos Package)
#if (ENABLE_STK_Mesh)
#    add_subdirectory(mesh_stk)
#endif(ENABLE_STK_Mesh)    

# MOAB Mesh framework
#if (ENABLE_MOAB_Mesh)
#    add_subdirectory(mesh_moab)
#endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
    add_subdirectory(mesh_mstk)
endif(ENABLE_MSTK_Mesh)    



# Factory class

add_subdirectory(mesh_factory)


if (BUILD_TESTS)

    # Include dirs for test
    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})

    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_SOURCE_DIR} EQUAL ${MESH_BINARY_DIR}))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${MESH_SOURCE_DIR}/test ${MESH_BINARY_DIR}/test) 
    endif()

    # Test link libraries
    set(test_link_libs geometry ${mesh_link_libs} simple_mesh mesh_factory ${UnitTest_LIBRARIES})

    # STK (Trilinos Package)
   if (ENABLE_STK_Mesh)
     list(APPEND test_link_libs stk_mesh)
   endif(ENABLE_STK_Mesh)    

    # MOAB Mesh framework
   if (ENABLE_MOAB_Mesh)
     list(APPEND test_link_libs moab_mesh)
   endif(ENABLE_MOAB_Mesh)    
    
    # MSTK Mesh framework
    if (ENABLE_MSTK_Mesh)
      list(APPEND test_link_libs mstk_mesh)
    endif(ENABLE_MSTK_Mesh)    

    # Test: mesh geometry

    add_Jali_test(mesh_geometry test_mesh_geometry
                    KIND unit
		    SOURCE test/Main.cc test/test_mesh_geometry.cc
		    LINK_LIBS ${test_link_libs})

    # Test mesh entities

    add_Jali_test(mesh_iterators test_entity_iterators
                  KIND unit
                  SOURCE test/Main.cc test/test_entity_iterators.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_entity_counts_serial mesh_entity_counts_serial
                  KIND unit
		  SOURCE test/Main.cc test/test_entity_counts.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_entity_counts_parallel mesh_entity_counts_parallel
                  KIND unit
		  NPROCS 4
		  SOURCE test/Main.cc test/test_entity_counts.cc
		  LINK_LIBS ${test_link_libs})

    # Test side entities

    add_Jali_test(side_tests test_sides
                  KIND unit
		  SOURCE test/Main.cc test/test_sides.cc
                  LINK_LIBS ${test_link_libs})

    # Test wedge entities

    add_Jali_test(wedge_tests test_wedges
                  KIND unit
		  SOURCE test/Main.cc test/test_wedges.cc
                  LINK_LIBS ${test_link_libs})

    # Test corner entities

    add_Jali_test(corner_tests test_corners
                  KIND unit
		  SOURCE test/Main.cc test/test_corners.cc
                  LINK_LIBS ${test_link_libs})

    # Test mesh tiles

    add_Jali_test(tile_tests test_one_tile
                  KIND unit
		  SOURCE test/Main.cc test/test_one_tile.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_tile_tests_serial test_meshtiles_serial
                  KIND unit
		  SOURCE test/Main.cc test/test_meshtiles.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_tile_tests_parallel test_meshtiles_parallel
                  KIND unit
		  NPROCS 4
		  SOURCE test/Main.cc test/test_meshtiles.cc
		  LINK_LIBS ${test_link_libs})

    # Test boundary ghosts

    add_Jali_test(mesh_boundary_ghost_tests_serial test_boundary_ghosts_serial
                  KIND unit
		  SOURCE test/Main.cc test/test_boundary_ghosts.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_boundary_ghost_tests_parallel test_boundary_ghosts_parallel
                  KIND unit
		  NPROCS 4
		  SOURCE test/Main.cc test/test_boundary_ghosts.cc
		  LINK_LIBS ${test_link_libs})


    # Test mesh sets

    add_Jali_test(mesh_set_tests_serial test_meshsets_serial
                  KIND unit
		  SOURCE test/Main.cc test/test_meshsets.cc
		  LINK_LIBS ${test_link_libs})

    add_Jali_test(mesh_set_tests_parallel test_meshsets_parallel
                  KIND unit
		  NPROCS 4
		  SOURCE test/Main.cc test/test_meshsets.cc
		  LINK_LIBS ${test_link_libs})
    

endif()
  
