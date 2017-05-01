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

# ############################################################################ #
#                                                                              #
#  DEFINE_EXTERNAL_PROJECT_ARGS(<PACK_NAME>                                    # 
#                                [ TARGET target-name ]                        #
#                                [ BUILD_IN_SOURCE ]                           #
#                                [ DEPENDS pack1 pack2 pack3 ]                 #
#                              )                                               #
#                                                                              #
#  A macro that defines common arguments for the AddExternalProject function   #
#  This macro provides an organized build structure.                           #
#                                                                              #
# ############################################################################ #
include(CMakeParseArguments)
include(SetMacros)
macro(DEFINE_EXTERNAL_PROJECT_ARGS prefix)
 
  # --- Parse the arguments
  set(_flags    "BUILD_IN_SOURCE")
  set(_oneValue "TARGET")
  set(_multiValue "DEPENDS")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})


  # --- Define the build target name
  if ( NOT PARSE_TARGET )
    string(TOLOWER "${prefix}" _target_name)
  else()
    set(_target_name ${PARSE_TARGET})
  endif()
  global_set(${prefix}_BUILD_TARGET ${_target_name})


  # --  Define the directories for download, build, install and timestamps

  # Will use lower case for directory names
  string(TOLOWER "${${prefix}_BUILD_TARGET}" target_lc) 


  set(${prefix}_prefix_dir ${SuperBuild_BINARY_DIR}/${target_lc})
  set(${prefix}_source_dir ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-${${prefix}_VERSION}-source)
  set(${prefix}_stamp_dir  ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-timestamps)
  set(${prefix}_tmp_dir    ${SuperBuild_BINARY_DIR}/${target_lc}/tmp)

  # Default is to build out of source, but some packages can not do that
  if ( NOT PARSE_BUILD_IN_SOURCE ) 
    set(${prefix}_build_dir  ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-${${prefix}_VERSION}-build)
  else()  
    set(${prefix}_build_dir "")
  endif()

  # Download from the web unless DISABLE_EXTERNAL_DOWNLOADS is TRUE
  if ( DISABLE_EXTERNAL_DOWNLOAD )
    set(${prefix}_URL ${TPL_DOWNLOAD_DIR}/${${prefix}_ARCHIVE_FILE})

    if ( NOT EXISTS "${${prefix}_URL}" )
      message(FATAL_ERROR "You have disabled external downloads (-DDISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE),"
	                  "however ${${prefix}_URL} does not exist")
    endif()

  else()
    set(${prefix}_URL ${${prefix}_URL_STRING}/${${prefix}_ARCHIVE_FILE})
  endif()

  # --- Set additional arguments

  # Log all steps, this keeps the STDOUT/STDERR tidy.
  set(${prefix}_logging_args
                          LOG_DOWNLOAD  1 
                          LOG_UPDATE    1 
                          LOG_CONFIGURE 1
                          LOG_BUILD     1
                          LOG_TEST      1
                          LOG_INSTALL   1)

  # Define the package dependencies 			
  set(${prefix}_PACKAGE_DEPENDS)
  foreach( _pack ${PARSE_DEPENDS})
    set(_pack_target "${${_pack}_BUILD_TARGET}")
    if ( NOT TARGET ${_pack_target} )
      message(FATAL_ERROR "Package ${prefix} requires ${_pack}, "
	                  "however the build configuration for ${_pack} has not been defined.")
    endif()
    list(APPEND ${prefix}_PACKAGE_DEPENDS "${_pack_target}")
  endforeach()

  # Set the build in source flag
  set(${prefix}_BUILD_IN_SOURCE ${PARSE_BUILD_IN_SOURCE})


endmacro(DEFINE_EXTERNAL_PROJECT_ARGS)
