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
#  jali
#    statemanager
#

# Jali module, include files found in JALI_MODULE_PATH
#include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   JALI_STATE_SOURCE_DIR
#   JALI_STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(JALI_STATE)

set(JALI_STATE_headers
  JaliState.h
  JaliStateVector.h
  )
list(TRANSFORM JALI_STATE_headers PREPEND "${JALI_STATE_SOURCE_DIR}/")

set(JALI_STATE_sources
  JaliState.cc
  JaliStateVector.cc
  )


#
# Library: 
#
add_library(jali_state ${JALI_STATE_sources})
set_target_properties(jali_state PROPERTIES PUBLIC_HEADER "${JALI_STATE_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_state ALIAS jali_state)

# We should change install_interface to include/error_handling but that requires
# changing all client codes
target_include_directories(jali_state PUBLIC
  $<BUILD_INTERFACE:${JALI_STATE_BINARY_DIR}>
  $<BUILD_INTERFACE:${JALI_STATE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

find_package(Boost QUIET COMPONENTS system filesystem program_options regex REQUIRED)
target_include_directories(jali_state PUBLIC ${Boost_INCLUDE_DIRS})

  
# Make the error handling and mesh targets a dependency of this target
target_link_libraries(jali_state PUBLIC jali_error_handling jali_mesh)

install(TARGETS jali_state
  EXPORT JaliTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

if (BUILD_TESTS)

  # Include dirs for test
#  target_link_libraries(jali_state PRIVATE mesh_factory)
#  target_link_libraries(jali_state PRIVATE mstk_mesh)
#  target_link_libraries(jali_state PRIVATE geometry)
#  target_link_libraries(jali_state PRIVATE ${UnitTest_LIBRARIES})

  # Need to copy files for the tests 
  if (NOT (${JALI_STATE_SOURCE_DIR} EQUAL
	${JALI_STATE_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${JALI_STATE_SOURCE_DIR}/test ${JALI_STATE_BINARY_DIR}/test) 
  endif()
    
  # Test state vectors
  
  set(test_src_files test/Main.cc test/test_jali_state_vector.cc)
  
  add_Jali_test(jali_state_vectors test_jali_state_vectors
    KIND unit
    SOURCE ${test_src_files}
    LINK_LIBS jali_state jali_mesh_factory ${UnitTest++_LIBRARIES})

  # Test state
  
  set(test_src_files test/Main.cc test/test_jali_state.cc)
  
  add_Jali_test(jali_state test_jali_state
    KIND unit
    SOURCE ${test_src_files}
    LINK_LIBS jali_state jali_mesh_factory ${UnitTest++_LIBRARIES})
endif()
  
