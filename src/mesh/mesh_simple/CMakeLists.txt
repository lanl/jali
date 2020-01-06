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
#    Mesh Base Class
#

# Jali module, include files found in JALI_MODULE_PATH
# include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SIMPLE_SOURCE_DIR
#   MESH_SIMPLE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_SIMPLE)

# Library: simple_mesh
set(MESH_SIMPLE_headers
  Mesh_simple.hh)
set(MESH_SIMPLE_sources
  Mesh_simple.cc)

add_library(jali_simple_mesh ${MESH_SIMPLE_headers} ${MESH_SIMPLE_sources})

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_simple_mesh ALIAS jali_simple_mesh)


target_include_directories(jali_simple_mesh PUBLIC
  $<BUILD_INTERFACE:${MESH_SIMPLE_BINARY_DIR}>
  $<BUILD_INTERFACE:${MESH_SIMPLE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(jali_simple_mesh PUBLIC jali_error_handling)
target_link_libraries(jali_simple_mesh PUBLIC jali_geometry)
target_link_libraries(jali_simple_mesh PUBLIC jali_mesh)

if (BUILD_TESTS)

  # Need to copy files for the tests 
  if (NOT (${MESH_SIMPLE_SOURCE_DIR} EQUAL ${MESH_SIMPLE_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${MESH_SIMPLE_SOURCE_DIR}/test ${MESH_SIMPLE_BINARY_DIR}/test) 
  endif()

  # Test: simple_mesh
  add_Jali_test(simple_mesh test_simple_mesh
    KIND unit
    SOURCE
    test/Main.cc 
    test/test_maps.cc 
    test/test_cell_numbering.cc 
    test/test_face_adj_cells.cc 
    test/test_node_adj_cells.cc 
    test/test_node_cell_faces.cc
    test/test_geometry.cc
    LINK_LIBS jali_simple_mesh ${UnitTest++_LIBRARIES})

endif()

