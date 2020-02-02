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
#    MSTK Mesh frameworks
#

# Jali module, include files found in JALI_MODULE_PATH
# include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MSTK_MESH_SOURCE_DIR
#   MSTK_MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MSTK_MESH)

set(MSTK_MESH_headers
  Mesh_MSTK.hh)
list(TRANSFORM MSTK_MESH_headers PREPEND "${MSTK_MESH_SOURCE_DIR}/")

set(MSTK_MESH_sources
  Mesh_MSTK.cc)

#
# Library: mstk_mesh
#

add_library(jali_mstk_mesh ${MSTK_MESH_sources})
set_target_properties(jali_mstk_mesh PROPERTIES PUBLIC_HEADER "${MSTK_MESH_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_mstk_mesh ALIAS jali_mstk_mesh)


# Find MSTK. MSTK needs to be built with ExodusII, METIS and ZOLTAN
find_package(MSTK QUIET REQUIRED)  # defines target mstk::mstk
set_property(TARGET ${MSTK_LIBRARIES} PROPERTY IMPORTED_GLOBAL TRUE)
if (MSTK_ROOT)
  set(MSTK_ROOT ${MSTK_ROOT} CACHE PATH "MSTK installation directory" FORCE)
else ()
  get_filename_component(MSTK_ROOT "${MSTK_INCLUDE_DIR}/.." ABSOLUTE)
  set(MSTK_ROOT ${MSTK_ROOT} CACHE PATH "MSTK installation directory" FORCE)
endif ()


target_link_libraries(jali_mstk_mesh PUBLIC ${MSTK_LIBRARIES})
target_link_libraries(jali_mstk_mesh PUBLIC jali_error_handling)
target_link_libraries(jali_mstk_mesh PUBLIC jali_mesh)

# We should change install_interface to include/jali/mstk_mesh but
# that requires changing all client codes
target_include_directories(jali_mstk_mesh PUBLIC
  $<BUILD_INTERFACE:${MSTK_MESH_BINARY_DIR}>
  $<BUILD_INTERFACE:${MSTK_MESH_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

install(TARGETS jali_mstk_mesh
  EXPORT JaliTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

if (BUILD_TESTS)
    
  # Need to copy files for the tests 
  if (NOT (${MSTK_MESH_SOURCE_DIR} EQUAL ${MSTK_MESH_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${MSTK_MESH_SOURCE_DIR}/test ${MSTK_MESH_BINARY_DIR}/test) 
  endif()
  

  # Test: mstk_mesh
  add_Jali_test(mstk_mesh test_mstk_mesh
    KIND int
    SOURCE
    test/Main.cc 
    test/test_hex_1x1x1.cc
    test/test_hex_2x2x1.cc
    test/test_hex_3x3x3.cc
    test/test_read_fractures.cc
    test/test_hex_gen_3x3x3.cc
    test/test_quad_gen_3x3.cc
    test/test_write_read_fields.cc
    test/test_block_partition.cc
    LINK_LIBS jali_mstk_mesh ${UnitTest++_LIBRARIES}) 
  
  # Test: mstk_mesh_parallel
  add_Jali_test(mstk_mesh_parallel test_parallel_mstk_mesh
    KIND unit
    NPROCS 4
    SOURCE
    test/Main.cc
    test/test_hex_3x3x3_4P.cc
    test/test_hex_3x3x3_par_read_4P.cc
    test/test_quad_gen_5x5_par.cc
    test/test_hex_gen_5x5x5_par.cc
    test/test_edges_4P.cc
    LINK_LIBS jali_mstk_mesh ${UnitTest++_LIBRARIES})
  
  # Test: mstk_mesh_parallel
  add_Jali_test(mstk_mesh_parallel8 test_parallel8_mstk_mesh
    KIND unit
    NPROCS 8
    SOURCE
    test/Main.cc
    test/test_quad_gen_5x5_par.cc
    test/test_hex_gen_5x5x5_par.cc
    LINK_LIBS jali_mstk_mesh ${UnitTest++_LIBRARIES})

endif()

