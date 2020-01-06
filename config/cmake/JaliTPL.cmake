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

# -*- mode: cmake -*-
# 
# Jali Third Party Library (TPL) Definitions
#

# Use <package_name>_ROOT variables in find_package (introduced in 3.12)
if (${CMAKE_VERSION} VERSION_GREATER 3.12)
  cmake_policy(SET CMP0074 NEW)
endif(${CMAKE_VERSION} VERSION_GREATER 3.12)

##############################################################################
# ------------------------ Required Libraries -------------------------------#
##############################################################################

##############################################################################
# MPI
##############################################################################
ENABLE_LANGUAGE(C CXX Fortran)  # Fortran for Exodus II fortran interface

find_package(MPI REQUIRED)
  
# Warn the user if MPI information is not found
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


##############################################################################
#-------------------------- Optional Libraries ------------------------------#
##############################################################################

##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################
option(ENABLE_UnitTest "Build Jali unit tests. Requires UnitTest++" ON)
if (ENABLE_UnitTest)
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
  else ()
    set(UnitTest++_DIR ${unittest_dir_save})
    find_package(UnitTest++ QUIET REQUIRED MODULE)  # fallback to module
  endif ()

  message(STATUS "Found UnitTest++: ${UnitTest++_LIBRARIES}")
endif()    


