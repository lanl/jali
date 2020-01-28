# Copyright (c) 2019, Triad National Security, LLC
# All rights reserved.
#
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
# Jali 
#       Root CMakeLists.txt file
#

cmake_minimum_required(VERSION 3.11)
set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND})  

# Define the project name
# This command will define
#     Jali_SOURCE_DIR
#     Jali_BINARY_DIR
project(Jali)

set(Jali_VERSION_MAJOR 1)
set(Jali_VERSION_MINOR 1)
set(Jali_VERSION_PATCH 0)


if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables
endif()


enable_language(C CXX)

# ENABLE C++11 support

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()



# Jali needs MPI

find_package(MPI REQUIRED)
  
if (MPI_FOUND)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "MPI C wrapper to use" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "MPI C++ wrapper to use" FORCE)
endif ()


# --- Jali uses MPI_EXEC* not MPIEXEC* variables. This allows the user to 
#     override the find package results.

# - MPI execute binary
if (MPIEXEC)
  set(MPI_EXEC ${MPIEXEC} CACHE STRING "Custom MPI Executable specified" FORCE)
else()
  set(MPI_EXEC ${MPIEXEC_EXECUTABLE} CACHE STRING "MPI Executable found by FindMPI" FORCE)
endif()


# INTERNAL - Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if (NOT MPI_EXEC_NUMPROCS_FLAG)
  if (MPIEXEC_NUMPROC_FLAG)
    set(MPI_EXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "Set MPI number of procs flag from FindMPI")
  else()
    set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
  endif()
endif()
  
# INTERNAL - Maximum number of processors for the test suite Some
# tests require too many processors and it increases the execution
# time considerably.
set(MPI_EXEC_MAX_NUMPROCS_DFLT 8)
if (NOT MPI_EXEC_MAX_NUMPROCS)
  include(ProcessorCount)
  ProcessorCount(proc_count)
  if (NOT proc_count EQUAL 0)
    math(EXPR MPI_EXEC_MAX_NUMPROCS "${proc_count} * 2") 
    message(STATUS "Detected ${proc_count} processors and will set maximum to ${MPI_EXEC_MAX_NUMPROCS}")
  else()
    set(MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DFLT})
  endif()
endif()  



# Useful variables pointing to directories in the source tree
set(Jali_SOURCE_TOOLS_DIR    "${Jali_SOURCE_DIR}/tools")
set(Jali_SOURCE_CONFIG_DIR   "${Jali_SOURCE_DIR}/config")
set(Jali_PYTHON_DIR          "${Jali_SOURCE_DIR}/tools/py_lib")

# Set the module search path so find_package and include commands
# can locate files in <root source tree>/tools/cmake
set(Jali_MODULE_PATH "${Jali_SOURCE_CONFIG_DIR}/cmake")
set(CMAKE_MODULE_PATH 
    ${Jali_MODULE_PATH}
    ${Jali_MODULE_PATH}/Modules)

# This will be the name of the overall target (no library named
# libJali will be built)

add_library(Jali INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::Jali ALIAS Jali)


# Mesh frameworks
option(ENABLE_MOAB_Mesh "Enable MOAB mesh framework" OFF)
option(ENABLE_STK_Mesh "Enable STK mesh framework" OFF)
option(ENABLE_MSTK_Mesh "Enable MSTK mesh framework" ON)
if (NOT ENABLE_MOAB_Mesh AND
    NOT ENABLE_STK_Mesh AND
    NOT ENABLE_MSTK_Mesh)
  message(STATUS "No parallel unstructured framework enabled?")
endif ()

# Testing
option(ENABLE_TESTS
  "Build Jali unit tests. Requires UnitTest++" ON)     # can be overridden

if (ENABLE_TESTS)
  # Enable testing, ctest needs this
  # all add_test commands are ignored unless this is called!
  enable_testing()
  
  set(BUILD_TESTS 1)

  # Some platforms require all binaries linking to MPI
  # only run through the MPIEXEC binary
  option(TESTS_REQUIRE_MPIEXEC "Run all tests with the MPIEXEC binary" FALSE)
  
  # Need this option if the PATH environment does not include '.'
  option(TESTS_REQUIRE_FULLPATH "Append full path to test binaries" TRUE)

  # For backward compatibility allow UnitTest_DIR for a bit more time
  if (NOT UnitTest++_DIR AND UnitTest_DIR)
    set(UnitTest++_DIR ${UnitTest_DIR})
  endif ()

  set(unittest_dir_save ${UnitTest++_DIR})
  find_package(UnitTest++ QUIET CONFIG PATHS ${UnitTest++_DIR})  # Try to discover thru cmake config file

  if (UnitTest++_FOUND)
    # UnitTest++ sets a weird path for it's includes
    set(UnitTest++_INCLUDE_DIRS ${UTPP_INCLUDE_DIRS}/UnitTest++)

    # Also it sets the target name as UnitTest++ but does not set the
    # LIBRARIES variable
    set(UnitTest++_LIBRARIES UnitTest++)

    # Finally it does not connect the include directories to the target
    set_target_properties(UnitTest++ PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${UnitTest++_INCLUDE_DIRS})
  else ()
    set(UnitTest++_DIR ${unittest_dir_save})
    find_package(UnitTest++ QUIET REQUIRED MODULE)  # fallback to module
  endif ()

  message(STATUS "Found UnitTest++: ${UnitTest++_LIBRARIES}")
endif()    


# Process source directories
# Recurse down and individual component libaries (error_handling,
# geometry, state, mesh)

add_subdirectory(src)

target_link_libraries(Jali INTERFACE jali_error_handling)
target_link_libraries(Jali INTERFACE jali_geometry)
target_link_libraries(Jali INTERFACE jali_mesh)
target_link_libraries(Jali INTERFACE jali_mesh_factory)
target_link_libraries(Jali INTERFACE jali_state)


# Example usage code for users

add_subdirectory(examples)


# Set up targets and export

set(Jali_LIBRARIES jali_error_handling jali_geometry jali_state jali_mesh jali_mesh_factory jali_simple_mesh)
if (ENABLE_MSTK_Mesh)
  list(APPEND Jali_LIBRARIES jali_mstk_mesh)
endif ()

install(TARGETS ${Jali_LIBRARIES}
  EXPORT JaliTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

install(EXPORT JaliTargets
  FILE JaliTargets.cmake
  NAMESPACE Jali::
  DESTINATION lib/cmake/Jali
  )
  
# Install the config file after filling in variables delimited by @..@
if (METIS_DIR AND NOT METIS_ROOT)
  set(METIS_ROOT ${METIS_DIR})
endif ()
if (Zoltan_DIR AND NOT Zoltan_ROOT)
  set(Zoltan_ROOT ${Zoltan_DIR})
endif ()
if (ExodusII_DIR AND NOT ExodusII_ROOT)
  set(ExodusII_ROOT ${ExodusII_DIR})
endif ()
configure_file(${PROJECT_SOURCE_DIR}/config/cmake/JaliConfig.cmake.in
  ${PROJECT_BINARY_DIR}/JaliConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/JaliConfig.cmake
  DESTINATION lib/cmake/Jali)

# Also install any module files used for finding the dependencies

install(FILES
  ${PROJECT_SOURCE_DIR}/config/cmake/Modules/FindMETIS.cmake
  ${PROJECT_SOURCE_DIR}/config/cmake/Modules/FindZoltan.cmake
  ${PROJECT_SOURCE_DIR}/config/cmake/Modules/FindExodusII.cmake
  DESTINATION lib/cmake/Jali/modules)

