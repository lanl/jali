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
# Functions for building and managing binaries.
#

include(CMakeParseArguments)
include(PrintVariable)
include(InstallManager)

#
# Usage:
#
# ADD_Jali_LIBRARY(<target name>
#                    SOURCE file1 file2  .....
#                    [HEADERS file1 file2 .....]
#                    [LINK_LIBS lib1 lib2 [link_opt1] ]
#                    [STATIC] [SHARED]
#                    [NO_INSTALL] [NO_INSTALL_HEADERS] )
# 
#
# Arguments:
#
#   target CMake target name defined for this library
#
#   SOURCE List of source files to compile 
#
#   HEADERS (Optional) List of heard files associated with this library
#
#   LINK_LIBS (Optional) Defines the list of link libraries required to build and link target
#
#   SHARED STATIC (Optional) Build a SHARED or STATIC library. If neither flag is set, then
#   CMake builds the library based on the BUILD_SHARED_LIBS setting. These flags can be used
#   to override BUILD_SHARED_LIBS.
#
#   NO_INSTALL (Optional) All libraries defined by this function will be added to the install
#   target. Use this option to not install the library.
#
#   NO_INSTALL_HEADERS (Optional) Any file found in the HEADERS variable will be added to the
#   install target. Use this option to not install the library.
#
# 
function(ADD_Jali_LIBRARY target)

  # --- Parse the input
  set(options STATIC SHARED NO_INSTALL NO_INSTALL_HEADERS)
  set(oneValueArgs "")
  set(multiValueArgs SOURCE HEADERS LINK_LIBS)
  cmake_parse_arguments(Jali_LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # --- Check the input
  if ( NOT Jali_LIB_SOURCE )
    message(FATAL "Can not define library ${target} without source files.")
  endif()

  if ( Jali_LIB_SHARED AND Jali_LIB_STATIC )
    message(FATAL "Library ${target} must either be STATIC or SHARED")
  endif()

  # --- Define the library
  if (Jali_LIB_SHARED)
    set(shared_flag SHARED)
  endif()
  if (Jali_LIB_STATIC)
    set(static_flag STATIC)
  endif()
  add_library(${target} ${shared_flag} ${static_flag} ${Jali_LIB_SOURCE})

  # --- Add link libraries
  if ( Jali_LIB_LINK_LIBS )
    target_link_libraries(${target} ${Jali_LIB_LINK_LIBS})
  endif()

  # --- Add target to the install target
  if ( NOT "${Jali_LIB_NO_INSTALL}" )
    add_install_library(${target})
  endif()

  # --- Add header files to the install target
  if ( NOT "${Jali_LIB_NO_INSTALL_HEADERS}" )
    add_install_include_file(${Jali_LIB_HEADERS})
  endif()


endfunction(ADD_Jali_LIBRARY)

