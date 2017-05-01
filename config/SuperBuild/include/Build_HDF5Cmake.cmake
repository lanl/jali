# Copyright (c) 2017, Los Alamos National Security, LLC
# All rights reserved.

# Copyright 2017. Los Alamos National Security, LLC. This software was
# produced under U.S. Government contract DE-AC52-06NA25396 for Los
# Alamos National Laboratory (LANL), which is operated by Los Alamos
# National Security, LLC for the U.S. Department of Energy. The
# U.S. Government has rights to use, reproduce, and distribute this
# software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
# LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
# derivative works, such modified software should be clearly marked, so
# as not to confuse it with the version available from LANL.
 
# Additionally, redistribution and use in source and binary forms, with
# or without modification, are permitted provided that the following
# conditions are met:

# 1.  Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2.  Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
# National Laboratory, LANL, the U.S. Government, nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
 
# THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
# ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#
# Build TPL:  HDF5 
#    
define_external_project(HDF5 TARGET hdf5 DEPENDS ZLIB )

# ############################################################################ #
# Create the patch shell script file
# ############################################################################ #
set(HDF5_sh_patch "${HDF5_prefix_dir}/hdf5-patch-step.sh")
configure_file(${SuperBuild_BUILD_FILES_DIR}/hdf5-patch-step.sh.in
               ${HDF5_sh_patch}
               @ONLY)

# ############################################################################ #
# Add external project
# ############################################################################ #
ExternalProject_Add(${HDF5_target}
    DEPENDS zlib
    ${HDF5_ep_directory_args}
    ${HDF5_url_args}
    PATCH_COMMAND sh ${HDF5_sh_patch}
    CMAKE_ARGS
             -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
             -D HDF5_ENABLE_PARALLEL:Bool=TRUE
             -D HDF5_ENABLE_Z_LIB_SUPPORT:bool=TRUE
             -D ZLIB_INCLUDE_DIR=${ZLIB_install_dir}/include
             -D ZLIB_LIBRARY=${ZLIB_install_dir}/lib/libz.a
             -D HDF5_BUILD_HL_LIB:bool=TRUE
             -D HDF5_BUILD_TOOLS:bool=TRUE
    ${HDF5_logging_opts}                  
)

