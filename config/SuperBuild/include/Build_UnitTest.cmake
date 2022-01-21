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
# Build TPL: UnitTest
# 
# --- Define all the directories and common external project flags
set(unittest_depend_projects ZLIB)

define_external_project_args(UnitTest
                             TARGET unittest
                             BUILD_IN_SOURCE)

# add version version to the autogenerated tpl_versions.h file
Jali_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX UnitTest
  VERSION ${UnitTest_VERSION_MAJOR} ${UnitTest_VERSION_MINOR} ${UnitTest_VERSION_PATCH})

# --- define the configuration parameters

set(Unittest_CMAKE_PACKAGE_ARGS "")

set(Unittest_CMAKE_TPL_ARGS)

# Pass the following MPI arguments to unittest
set(MPI_CMAKE_ARGS DIR EXEC EXEC_NUMPROCS_FLAG EXE_MAX_NUMPROCS C_COMPILER)
foreach (var ${MPI_CMAKE_ARGS} )
  set(mpi_var "MPI_${var}")
  if ( ${mpi_var} )
    list(APPEND Unittest_CMAKE_TPL_ARGS "-D${mpi_var}:STRING=${${mpi_var}}")
  endif()
endforeach() 

# build type
if ( CMAKE_BUILD_TYPE )
  list(APPEND Unittest_CMAKE_EXTRA_ARGS
              "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
endif()

# - Final language ARGS
set(Unittest_CMAKE_LANG_ARGS
                   ${Jali_CMAKE_C_COMPILER_ARGS}
		           -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                   ${Jali_CMAKE_CXX_COMPILER_ARGS}
		           -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                   ${Jali_CMAKE_Fortran_COMPILER_ARGS}
                   -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER})


#  - Add CMake configuration file
if(Unittest_Build_Config_File)
    list(APPEND Unittest_Config_File_ARGS
        "-C${Unittest_Build_Config_File}")
    message(STATUS "Will add ${Unittest_Build_Config_File} to the Unittest configure")    
    message(DEBUG "Unittest_CMAKE_EXTRA_ARGS = ${Unittest_CMAKE_EXTRA_ARGS}")
endif()    

#  - Final CMake Arguments 
set(Unittest_CMAKE_ARGS 
   ${Unittest_CMAKE_PACKAGE_ARGS}
   ${Unittest_CMAKE_TPL_ARGS}
   ${Unittest_CMAKE_EXTRA_ARGS}
   ${Unittest_CMAKE_LANG_ARGS})


# --- Add external project build and tie to the ZLIB build target
ExternalProject_add(${UnitTest_BUILD_TARGET}
                    DEPENDS   ${UnitTest_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${UnitTest_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${UnitTest_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                  # Download directory
                    URL          ${UnitTest_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${UnitTest_MD5_SUM}                  # md5sum of the archive file
		    # -- Configure
		    SOURCE_DIR          ${UnitTest_source_dir}        # Defining forces CMake to mkdir SOURCE_DIR
		    CMAKE_ARGS          ${Unittest_Config_File_ARGS}
                    CMAKE_CACHE_ARGS    ${Unittest_CMAKE_ARGS}
                                        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		    # -- Build
		    BUILD_COMMAND       $(MAKE)                       # Run make in build directory $(MAKE) enables parallel build
		    BINARY_DIR          ${UnitTest_build_dir}         # Define the build directory
		    BUILD_IN_SOURCE     ${UnitTest_BUILD_IN_SOURCE}   # Flag in/out source build
                    # -- Install
                    INSTALL_DIR         ${CMAKE_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${UnitTest_logging_args})

