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
# Build TPL: SEACAS 
#    
# --- Define all the directories and common external project flags

# SEACAS does not call MPI directly, however HDF5 requires
# MPI and to resolve links we need MPI compile wrappers.
define_external_project_args(SEACAS
                             TARGET seacas
                             DEPENDS HDF5 NetCDF)


# add version version to the autogenerated tpl_versions.h file
Jali_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX SEACAS
  VERSION ${SEACAS_VERSION_MAJOR} ${SEACAS_VERSION_MINOR} ${SEACAS_VERSION_PATCH})
  
# --- Define the configure parameters

# Compile flags
set(seacas_cflags_list -I${CMAKE_INSTALL_PREFIX}/include ${Jali_COMMON_CFLAGS})
build_whitespace_string(seacas_cflags ${seacas_cflags_list})

set(seacas_cxxflags_list -I${CMAKE_INSTALL_PREFIX}/include ${Jali_COMMON_CXXFLAGS})
build_whitespace_string(seacas_cflags ${seacas_cxxflags_list})

set(seacas_fcflags_list -I${CMAKE_INSTALL_PREFIX}/include ${Jali_COMMON_FCFLAGS})
build_whitespace_string(seacas_fcflags ${seacas_fcflags_list})

set(seacas_lflags_list)
build_whitespace_string(seacas_lflags ${seacas_lflags_list})

# Build the NetCDF libraries string
include(BuildLibraryName)
build_library_name(netcdf seacas_netcdf_library STATIC APPEND_PATH ${CMAKE_INSTALL_PREFIX}/lib)
build_library_name(hdf5_hl seacas_hdf5_hl_library STATIC APPEND_PATH ${CMAKE_INSTALL_PREFIX}/lib)
build_library_name(hdf5 seacas_hdf5_library STATIC APPEND_PATH ${CMAKE_INSTALL_PREFIX}/lib)
build_library_name(z seacas_z_library STATIC APPEND_PATH ${CMAKE_INSTALL_PREFIX}/lib)
set(seacas_netcdf_libraries
       ${seacas_netcdf_library}
       ${seacas_hdf5_hl_library}
       ${seacas_hdf5_library} -ldl
       ${seacas_z_library})
if ((NOT MPI_WRAPPERS_IN_USE) AND (MPI_C_LIBRARIES))
  list(APPEND seacas_netcdf_libraries ${MPI_C_LIBRARIES})
endif()

#
# --- Define the SEACAS patch step - mainly for nem_slice to be able
# --- to handle columns
#
set(ENABLE_SEACAS_Patch ON)
if (ENABLE_SEACAS_Patch)
  set(SEACAS_patch_file seacas-nemslice.patch)
  configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/seacas-patch-step.sh.in
                 ${SEACAS_prefix_dir}/seacas-patch-step.sh
                 @ONLY)
  set(SEACAS_PATCH_COMMAND bash ${SEACAS_prefix_dir}/seacas-patch-step.sh)
  message(STATUS "Applying SEACAS patches")
else (ENABLE_SEACAS_Patch)
  set(SEACAS_PATCH_COMMAND)
  message(STATUS "Patch NOT APPLIED for SEACAS")
endif (ENABLE_SEACAS_Patch)

# --- Configure the package
set(SEACAS_CMAKE_CACHE_ARGS
                    -DCMAKE_INSTALL_PREFIX:FILEPATH=<INSTALL_DIR>
                    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                    -DCMAKE_EXE_LINKER_FLAGS:STRING=${seacas_lflags}
                    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
                    -DSEACASProj_ENABLE_ALL_PACKAGES:BOOL=FALSE
                    -DSEACASProj_ENABLE_SEACASExodus:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASNemslice:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASNemspread:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASExodiff:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASExotxt:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASExoformat:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_SEACASDecomp:STRING=:BOOL=TRUE
                    -DSEACASProj_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=FALSE
                    -DSEACASProj_ENABLE_SECONDARY_TESTED_CODE:BOOL=FALSE
                    -DSEACASProj_ENABLE_TESTS:BOOL=FALSE
                    -DSEACASProj_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON 
                    -DSEACASProj_HIDE_DEPRECATED_CODE:STRING="NO"
                    -DTPL_ENABLE_Netcdf:BOOL=TRUE
                    -DTPL_Netcdf_LIBRARIES:STRING=${seacas_netcdf_libraries}
                    -DNetcdf_INCLUDE_DIRS:STRING=${CMAKE_INSTALL_PREFIX}/include
                    -DTPL_Netcdf_PARALLEL:BOOL=TRUE
                    -DTPL_ENABLE_Matio:BOOL=FALSE
                    -DTPL_ENABLE_X11:BOOL=FALSE
                    -DTPL_ENABLE_CGNS:BOOL=FALSE
                    -DTPL_ENABLE_MPI:BOOL=ON
                    -DTPL_ENABLE_Pamgen:BOOL=FALSE
                    -DTPL_ENABLE_Pthread:BOOL=FALSE
                    -DSEACASExodus_ENABLE_THREADSAFE:BOOL=OFF
                    -DSEACASIoss_ENABLE_THREADSAFE:BOOL=OFF
                    -DCMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_PREFIX}/SEACAS/lib
                    -DCMAKE_INSTALL_NAME_DIR:PATH=${CMAKE_INSTALL_PREFIX}/SEACAS/lib)

# --- Add external project build and tie to the SEACAS build target
ExternalProject_Add(${SEACAS_BUILD_TARGET}
                    DEPENDS   ${SEACAS_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${SEACAS_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${SEACAS_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                # Download directory
                    URL          ${SEACAS_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${SEACAS_MD5_SUM}                  # md5sum of the archive file
                    # -- Patch
                    PATCH_COMMAND ${SEACAS_PATCH_COMMAND}
                    # -- Configure
                    SOURCE_DIR       ${SEACAS_source_dir}           # Source directory
                    CMAKE_CACHE_ARGS ${SEACAS_CMAKE_CACHE_ARGS}
                                     -DCMAKE_C_FLAGS:STRING=${Jali_COMMON_CFLAGS}  # Ensure uniform build
                                     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                                     -DCMAKE_CXX_FLAGS:STRING=${Jali_COMMON_CXXFLAGS}
                                     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                                     -DCMAKE_Fortran_FLAGS:STRING=${Jali_COMMON_FCFLAGS}
                                     -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
                    # -- Build
                    BINARY_DIR        ${SEACAS_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${SEACAS_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}/SEACAS   # Install directory, NOT in the usual place!
                    # -- Output control
                    ${SEACAS_logging_args})

# --- Useful variables for other packages that depend on SEACAS
set(SEACAS_DIR ${CMAKE_INSTALL_PREFIX}/SEACAS)
