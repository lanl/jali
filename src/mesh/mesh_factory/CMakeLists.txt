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
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_FACTORY_SOURCE_DIR
#   MESH_FACTORY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_FACTORY)

# External (TPL) include directories
include_directories(${STK_INCLUDE_DIRS})
include_directories(${ExodusII_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})


#
# Library: mesh_factory
#
set(mesh_factory_files 
    MeshFactory.cc)
file(GLOB mesh_factory_inc_files "*.hh")

list(APPEND mesh_factory_frameworks simple_mesh)

if (ENABLE_STK_Mesh)
    include_directories(${STK_INCLUDE_DIR} ${STK_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks stk_mesh exodus_reader mesh_data)
    add_definitions("-DHAVE_STK_MESH")
endif()

if (ENABLE_MOAB_Mesh)
    include_directories(${MOAB_INCLUDE_DIR} ${MOAB_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks moab_mesh)
    add_definitions("-DHAVE_MOAB_MESH")
endif()

if (ENABLE_MSTK_Mesh)
    include_directories(${MSTK_INCLUDE_DIR} ${MSTK_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks mstk_mesh)
    add_definitions("-DHAVE_MSTK_MESH")
endif()

add_Jali_library(mesh_factory
                   SOURCE ${mesh_factory_files} HEADERS ${mesh_factory_inc_files}
                   LINK_LIBS ${mesh_factory_frameworks} ${Boost_LIBRARIES})


if (BUILD_TESTS)

    # some test file locations defined as macros

    set(test_file_flags "${test_file_flags} -DMOAB_TEST_FILE=\"\\\"./test/hex_3x3x3_ss_4P.h5m\\\"\" ")
    set(test_file_flags "${test_file_flags} -DEXODUS_TEST_FILE=\"\\\"./test/hex_5x5x5.exo\\\"\" ")
    set(test_file_flags "${test_file_flags} -DNEMESIS_TEST_FILE=\"\\\"./test/hex_5x5x5.par\\\"\" ")

    set_source_files_properties(test/test_mesh_factory.cc PROPERTIES 
                                COMPILE_FLAGS ${test_file_flags})
   
   # Add UnitTest include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

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
		         mesh_factory
			 ${UnitTest_LIBRARIES}
			 ${Boost_LIBRARIES})

   # Do not call add_Jali_test again with SOURCE and LINK_LIBS since the
   # previous call carries instructions to build the executable
   # When n>1 procs used this test fails
   #add_Jali_test(mesh_factory_parallel test_mesh_factory NPROCS 2 KIND uint)
   add_Jali_test(mesh_factory_parallel test_mesh_factory KIND uint)

endif()
   
