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


# file: CMakeLists.txt
#
# Jali
#   Mesh Factory
# -------------------------------------------------------------

# Jali module, include files found in JALI_MODULE_PATH
include(TestManager)
# include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_FACTORY_SOURCE_DIR
#   MESH_FACTORY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_FACTORY)



#
# Library: mesh_factory
#
set(MESH_FACTORY_headers
  MeshFactory.hh)
list(TRANSFORM MESH_FACTORY_headers PREPEND "${MESH_FACTORY_SOURCE_DIR}/")

set(MESH_FACTORY_sources 
  MeshFactory.cc)

add_library(jali_mesh_factory ${MESH_FACTORY_sources})
set_target_properties(jali_mesh_factory PROPERTIES PUBLIC_HEADER "${MESH_FACTORY_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_mesh_factory ALIAS jali_mesh_factory)


target_include_directories(jali_mesh_factory PUBLIC
  $<BUILD_INTERFACE:${MESH_FACTORY_BINARY_DIR}>
  $<BUILD_INTERFACE:${MESH_FACTORY_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(jali_mesh_factory PUBLIC jali_error_handling)
target_link_libraries(jali_mesh_factory PUBLIC jali_mesh)


install(TARGETS jali_mesh_factory
  EXPORT JaliTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )


if (BUILD_TESTS)

  # some test file locations defined as macros
  
  set(test_file_flags "${test_file_flags} -DMOAB_TEST_FILE=\"\\\"./test/hex_3x3x3_ss_4P.h5m\\\"\" ")
  set(test_file_flags "${test_file_flags} -DEXODUS_TEST_FILE=\"\\\"./test/hex_5x5x5.exo\\\"\" ")
  set(test_file_flags "${test_file_flags} -DNEMESIS_TEST_FILE=\"\\\"./test/hex_5x5x5.par\\\"\" ")
  
  set_source_files_properties(test/test_mesh_factory.cc PROPERTIES 
                              COMPILE_FLAGS ${test_file_flags})
  
  # Add UnitTest include directories
  
  # Need to copy files for the tests 
  if (NOT (${MESH_FACTORY_SOURCE_DIR} EQUAL ${MESH_FACTORY_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${MESH_FACTORY_SOURCE_DIR}/test ${MESH_FACTORY_BINARY_DIR}/test) 
  endif()

   # Test: test_mesh_factory
   add_Jali_test(mesh_factory test_mesh_factory
     KIND unit
     SOURCE
     test/Main.cc
     test/test_mesh_factory.cc
     LINK_LIBS
     jali_mesh_factory ${UnitTest++_LIBRARIES})

   add_Jali_test(mesh_factory_parallel test_mesh_factory_parallel
     KIND unit
     SOURCE
     test/Main.cc
     test/test_mesh_factory.cc
     LINK_LIBS
     jali_mesh_factory ${UnitTest++_LIBRARIES}
     NPROCS 4)

endif()
   
