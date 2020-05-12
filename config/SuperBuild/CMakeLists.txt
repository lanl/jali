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
# Jali SuperBuild
#       Jali TPL Build Targets
#
# Set minimum required version of CMake
# NB: Must go before project() call.
#
cmake_minimum_required(VERSION 3.11)

#
# Language standard
# NB: Seems like a good idea to have this before project() call.
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ############################################################################ #
# Define project name (SuperBuild_SOURCE_DIR and SuperBuild_BINARY_DIR)
# ############################################################################ #
message(STATUS "\n------------------------------------------------")
message(STATUS "\n-- CMake: Configuring Jali TPLs build/install.\n--")
message(STATUS "---------------------------------------------")
project(SuperBuild C CXX Fortran)


# ############################################################################ #
# TPL versions
# ############################################################################ #
include(${SuperBuild_SOURCE_DIR}/TPLVersions.cmake)

# prepare the version numbers file
set(TPL_VERSIONS_INCLUDE_FILE ${CMAKE_BINARY_DIR}/tpl_versions.h)
file(WRITE ${TPL_VERSIONS_INCLUDE_FILE} "")
install(FILES ${CMAKE_BINARY_DIR}/tpl_versions.h DESTINATION include)
Jali_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX JALI
  VERSION ${JALI_TPLS_VERSION_MAJOR} ${JALI_TPLS_VERSION_MINOR} ${JALI_TPLS_VERSION_PATCH})

# ############################################################################ #
# Define Jali directories 
# ############################################################################ #
file(TO_CMAKE_PATH ${SuperBuild_SOURCE_DIR}/../.. Jali_SOURCE_DIR)

# ############################################################################ #
# Update the CMake module path
# ############################################################################ #
set(Jali_MODULE_PATH
  ${Jali_SOURCE_DIR}/config/cmake
  ${Jali_SOURCE_DIR}/config/cmake/Utils
  ${Jali_SOURCE_DIR}/config/cmake/Modules)

set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} 
  ${SuperBuild_SOURCE_DIR}/cmake
  ${Jali_MODULE_PATH})

# ############################################################################ #
# Location of the Build_* files and template files
# ############################################################################ #
set(SuperBuild_BUILD_FILES_DIR    "${SuperBuild_SOURCE_DIR}/include")
set(SuperBuild_TEMPLATE_FILES_DIR "${SuperBuild_SOURCE_DIR}/templates")


# ############################################################################ #
# Build options
# ############################################################################ #
set(SuperBuild_BUILD_TYPE_DFLT "Release")
if (NOT CMAKE_BUILD_TYPE) 
  set(CMAKE_BUILD_TYPE ${SuperBuild_BUILD_TYPE_DFLT})
endif()
message(STATUS "Setting build type -- ${CMAKE_BUILD_TYPE}")

# Download directory TPL_DOWNLOAD_DIR
set(TPL_DOWNLOAD_DIR_DFLT ${SuperBuild_BINARY_DIR}/Downloads)
if (NOT TPL_DOWNLOAD_DIR)
  set(TPL_DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR_DFLT})
endif()
message(STATUS "Download TPL files in ${TPL_DOWNLOAD_DIR}") 

# Disable external web downloads
option(DISABLE_EXTERNAL_DOWNLOAD "Flag to disable external web downloads" FALSE)
if (DISABLE_EXTERNAL_DOWNLOAD)
  message(STATUS "Will not download files from the web, instead will search "
                 "for archive files in ${TPL_DOWNLOAD_DIR}")
else()
  include(CheckDownload)
  check_download(TEST_FILE MD5SUMS TEST_URL ${JALI_TPLS_DOWNLOAD_URL})
endif()


# Jali unstructured mesh API
option(ENABLE_Unstructured "Enable Jali unstructured mesh capability" TRUE)

# Jali unstructured mesh API and frameworks
option(ENABLE_MSTK_Mesh "Enable the Unstructured MSTK Mesh Toolkit" TRUE)
option(ENABLE_STK_Mesh "Enable the Unstructured STK Mesh Toolkit" FALSE)
option(ENABLE_MOAB_Mesh "Enable the Unstructured MOAB Mesh Toolkit" FALSE)

#option(ENABLE_OpenMP "Add OpenMP directives to the build" TRUE)
set (ENABLE_OpenMP FALSE)
if ( ENABLE_OpenMP )
  find_package(OpenMP)
  find_package(OpenMP_Fortran)
  if ( OPENMP_FOUND )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OpenMP_LINK_FLAGS}")
  else()  
    message(SEND_ERROR "Failed to determine the OpenMP compiler flags. Will deactivate.")
    set(ENABLE_OpenMP FALSE CACHE BOOL "Enable OpenMP compile directives" FORCE)
  endif()
endif()

# If unstructured is not enabled, there is no reason to have these mesh frameworks
if ( NOT ENABLE_Unstructured )
  set (ENABLE_STK_Mesh FALSE)
  set (ENABLE_MSTK_Mesh FALSE)
  set (ENABLE_MOAB_Mesh FALSE)
endif()

# Jali Test Suite
option(ENABLE_TESTS "Build the Jali Test Suite" TRUE)


# ############################################################################ #
# Compiler checks
# ############################################################################ #
include(FortranCInterface)
FortranCInterface_VERIFY()

# ############################################################################ #
# MPI Configuration 
# ############################################################################ #


find_package(MPI REQUIRED)
  
# Warn the user if MPI information is not found
if (NOT MPI_C_FOUND)
  message(WARNING "Failed to locate MPI C include and library files")
endif()
if (NOT MPI_CXX_FOUND)
  message(WARNING "Failed to locate MPI C++ include and library files")
endif()

# We have to compile Jali and Jali TPLs with MPI - so forcibly set the compilers to the MPI wrappers 
set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "MPI C wrapper to use" FORCE)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "MPI C++ wrapper to use" FORCE)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE FILEPATH "MPI Fortran wrapper to use" FORCE)
  
# The FindMPI package uses FindPackageHandleStandardArgs
# to MPI_<lang>_FOUND. For consistency, set MPI_FORTRAN_FOUND
# to MPI_Fortran_FOUND
set(MPI_Fortran_FOUND ${MPI_FORTRAN_FOUND})

if (NOT MPI_Fortran_FOUND)
  message(WARNING "Failed to locate MPI Fortran include and library files")
endif()

# - Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if (NOT MPI_EXEC_NUMPROCS_FLAG)
  if (MPIEXEC_NUMPROC_FLAG)
    set(MPI_EXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "Set MPI number of procs flag from FindMPI")
  else()
    set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
  endif()
endif()
  
# - Maximum number of processors. This is a limit for the test suite
#   Some tests require too many processors and it increases the execution time
#   considerably. 
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

# - Set the pre and post flags
#   Usage:
#   ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} PROCS ${MPI_EXEC_PREFLAGS} EXECUTABLE ${MPI_EXEC_POSTFLAGS}
if (NOT MPI_EXEC_PREFLAGS)
  if (MPIEXEC_PREFLAGS)
    set(MPI_EXEC_PREFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute pre flags")
  endif()
endif()

if (NOT MPI_EXEC_POSTFLAGS)
  if (MPIEXEC_POSTFLAGS)
    set(MPI_EXEC_POSTFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute post flags")
  endif()
endif()


# ############################################################################ #
# Set common build compiler flags, build types and directories
# ############################################################################ #

# ############################################################################ #
# Shared Libraries Definitions
# ############################################################################ #
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 

set(JALI_CMAKE_CACHE_ARGS  "-DCMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_RPATH}")
list(APPEND JALI_CMAKE_CACHE_ARGS  "-DCMAKE_INSTALL_NAME_DIR:PATH=${CMAKE_INSTALL_RPATH}")

if (BUILD_SHARED_LIBS)
  set(_shared_flags "-fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_shared_flags}")

  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_shared_flags}")

  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_shared_flags}")
endif()

# ############################################################################ #
# Static Libraries Definitions
# ############################################################################ #
set(BUILD_STATIC_LIBS TRUE)
if (BUILD_SHARED_LIBS)
  set(BUILD_STATIC_LIBS FALSE)
endif()


# CMake compiler settings for any package built with CMake

# C language flags
set(Jali_CMAKE_C_COMPILER_ARGS
     -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
     -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
     -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
     -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
     -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO})

# C++ language flags 
set(Jali_CMAKE_CXX_COMPILER_ARGS
     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
     -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
     -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
     -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
     -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

# Fortran language flags    
set(Jali_CMAKE_Fortran_COMPILER_ARGS
     -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
     -DCMAKE_Fortran_FLAGS_DEBUG:STRING=${CMAKE_Fortran_FLAGS_DEBUG}
     -DCMAKE_Fortran_FLAGS_MINSIZEREL:STRING=${CMAKE_Fortran_FLAGS_MINSIZEREL}
     -DCMAKE_Fortran_FLAGS_RELEASE:STRING=${CMAKE_Fortran_FLAGS_RELEASE}
     -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})

# Link flags
set(Jali_CMAKE_LINKER_ARGS
     -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
     -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS})

# ############################################################################ #
# Common compiler flags to ensure a uniform build in projects that do not use CMake
# ############################################################################ #
include(BuildWhitespaceString)
string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_uc)
set(build_c_flags ${CMAKE_C_FLAGS_${build_type_uc}})
set(build_cxx_flags ${CMAKE_CXX_FLAGS_${build_type_uc}})
set(build_fortran_flags ${CMAKE_Fortran_FLAGS_${build_type_uc}})
message(STATUS "Common compile flags for build type ${CMAKE_BUILD_TYPE}")
build_whitespace_string(Jali_COMMON_CFLAGS ${CMAKE_C_FLAGS} ${build_c_flags}) 
build_whitespace_string(Jali_COMMON_CXXFLAGS ${CMAKE_CXX_FLAGS} ${build_cxx_flags}) 
build_whitespace_string(Jali_COMMON_FCFLAGS ${CMAKE_Fortran_FLAGS} ${build_fortran_flags}) 

# Common link flags for projects that do not use CMake
set(Jali_COMMON_LDFLAGS)
if (CMAKE_EXE_LINK_FLAGS)
  build_whitespace_string(Jali_COMMON_LDFLAGS ${CMAKE_EXE_LINK_FLAGS})
endif()

# GNU configure static/shared
if(BUILD_SHARED_LIBS)
  set(Jali_SHARED_SWITCH "--enable-shared")
else()
  set(Jali_SHARED_SWITCH "--disable-shared")
endif()

message(STATUS "\tcommon C flags\t\t${Jali_COMMON_CFLAGS}")
message(STATUS "\tcommon C++ flags\t${Jali_COMMON_CXXFLAGS}")
message(STATUS "\tcommon Fortran flags\t${Jali_COMMON_FCFLAGS}")
message(STATUS "\tcommon Linker flags\t${Jali_COMMON_LDFLAGS}")


# ############################################################################ #
# LAPACK/BLAS Configuration 
# ############################################################################ #

# Search for BLAS (Basic Linear Algebra Subprograms) http://www.netlib.org/blas
find_package(BLAS REQUIRED)

# Search for LAPACK (Linear Algebra PACKage) http://www/netlib.org/lapack
find_package(LAPACK REQUIRED)

# ############################################################################ #
# Begin TPL builds 
# ############################################################################ #

# --- Add .NOTPARALLEL target in order to force the top level make to 
#     execute serially

add_custom_target(.NOTPARALLEL)

# --- Initalize and set the TPL build files

# Create configure, build, install and test targets for each TPL
set_property(DIRECTORY PROPERTY
             EP_STEP_TARGETS download patch configure build install test)

# Include the TPL version information
include(${SuperBuild_SOURCE_DIR}/TPLVersions.cmake)

# Include the CMake module ExternalProject
include(ExternalProject)

# Include the macro that defines uniform build, install and Add_ExternalProject args
include(DefineExternalProjectArgs)

# Set of CMake configure arguments to pass to Jali once the TPLs are built
set(Jali_TPL_CMAKE_ARGS)

# --- Begin the build definitions

# MPI
if ( BUILD_MPI )
  include(${SuperBuild_BUILD_FILES_DIR}/Build_OpenMPI.cmake)
  set(MPI_PROJECT OpenMPI)
endif()

# ZLIB
include(${SuperBuild_BUILD_FILES_DIR}/Build_zlib.cmake) 

# UnitTest
append_set(Jali_TPL_CMAKE_ARGS
           -DENABLE_TEST:BOOL=${ENABLE_TESTS})
if (ENABLE_TESTS)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_UnitTest.cmake)
  append_set(Jali_TPL_CMAKE_ARGS 
             -DUnitTest_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})
endif()

# Boost
include(${SuperBuild_BUILD_FILES_DIR}/Build_Boost.cmake) 
append_set(Jali_TPL_CMAKE_ARGS
           -DBOOST_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX}
           -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE)

# HDF5  
include(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5.cmake) 
append_set(Jali_TPL_CMAKE_ARGS
           -DHDF5_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})

# NetCDF
include(${SuperBuild_BUILD_FILES_DIR}/Build_NetCDF.cmake) 
append_set(Jali_TPL_CMAKE_ARGS
           -DNetCDF_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})

# NetCDF Fortran
#include(${SuperBuild_BUILD_FILES_DIR}/Build_NetCDF_Fortran.cmake)
#append_set(Jali_TPL_CMAKE_ARGS
#           -DNetCDF_FORTRAN_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})

# ExodusII
include(${SuperBuild_BUILD_FILES_DIR}/Build_ExodusII.cmake) 
append_set(Jali_TPL_CMAKE_ARGS
           -DExodusII_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})

# MOAB
if ( ENABLE_MOAB_Mesh )
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MOAB.cmake)
  append_set(Jali_TPL_CMAKE_ARGS
             -DENABLE_MOAB_Mesh:BOOL=${ENABLE_MOAB_Mesh}
	     -DMOAB_ROOT:PATH=${CMAKE_INSTALL_PREFIX})
endif()

include(${SuperBuild_BUILD_FILES_DIR}/Build_METIS.cmake) 
append_set(Jali_TPL_CMAKE_ARGS 
           -DMETIS_ROOT:FILEPATH=${CMAKE_INSTALL_PREFIX})

# Trilinos  
include(${SuperBuild_BUILD_FILES_DIR}/Build_Trilinos.cmake)
append_set(Jali_TPL_CMAKE_ARGS
           -DENABLE_STK_Mesh:BOOL=${ENABLE_STK_Mesh}
           -DTrilinos_INSTALL_PREFIX:PATH=${Trilinos_INSTALL_PREFIX}
	   -DTrilinos_ROOT:FILEPATH=${Trilinos_INSTALL_PREFIX})

# SEACAS
#include(${SuperBuild_BUILD_FILES_DIR}/Build_SEACAS.cmake)

# MSTK
if ( ENABLE_MSTK_Mesh )
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MSTK.cmake) 
  append_set(Jali_TPL_CMAKE_ARGS
            -DENABLE_MSTK_Mesh:BOOL=${ENABLE_MSTK_Mesh}
            -DMSTK_ROOT:PATH=${CMAKE_INSTALL_PREFIX})
endif()

# ############################################################################ #
#  Final steps
# ############################################################################ #

#  --- Create a CMake cache scipt
set(SuperBuild_Jali_CACHE_SCRIPT Jali-tpl-config.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/${SuperBuild_Jali_CACHE_SCRIPT}.in
               ${SuperBuild_BINARY_DIR}/${SuperBuild_Jali_CACHE_SCRIPT}
               @ONLY)
install( FILES "${SuperBuild_BINARY_DIR}/${SuperBuild_Jali_CACHE_SCRIPT}"             
         DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake)

#  --- Build the docs
option(BUILD_DOCS "Build the SuperBuild documentation" FALSE)
if ( BUILD_DOCS )
  add_subdirectory(doc)
endif()

