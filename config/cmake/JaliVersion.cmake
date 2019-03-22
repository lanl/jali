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
# Functions for building and managing binaries.
#

include(CMakeParseArguments)
include(PrintVariable)
include(InstallManager)

#
# Usage:
#
# ADD_Jali_EXECUTABLE(<target name>
#                       SOURCE file1 file2 ....
#                       LINK_LIBS lib1 lib2 .....
#                       OUPUT_NAME name
#                       OUTPUT_DIRECTORY dir
#                       NO_INSTALL)
# 
#
# Arguments:
#
#   target CMake target name defined for this binary
#
#   SOURCE List of source files to compile 
#
#   LINK_LIBS Defines the list of link libraries required to build and link target
#
#   OUTPUT_NAME (Optional) Set the output name to OUTPUT_NAME, will use the CMake defaults
#   if this is not set.
#
#   OUTPUT_DIRECTORY (Optional) Set the output directory, default is CMAKE_CURRENT_BINARY_DIR
#
#   NO_INSTALL (Optional) All binaries defined by this function will be added to the install
#   target. Use this option to not install the binary.
#
# 
function(ADD_Jali_EXECUTABLE target)

  # --- Parse the input
  set(options NO_INSTALL)
  set(oneValueArgs OUTPUT_NAME OUTPUT_DIRECTORY)
  set(multiValueArgs SOURCE LINK_LIBS)
  cmake_parse_arguments(Jali_BIN "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # --- Check the input
  if ( NOT Jali_BIN_SOURCE )
    message(FATAL "Can not define binary ${target} without source files.")
  endif()

  # --- Define the executable
  add_executable(${target} ${Jali_BIN_SOURCE})

  # --- Add link libraries
  if ( Jali_BIN_LINK_LIBS )
    target_link_libraries(${target} ${Jali_BIN_LINK_LIBS})
  endif()

  # --- Add target to the install target
  if ( NOT "${Jali_BIN_NO_INSTALL}" )
    add_install_binary(${target})
  endif()

  # --- Change the output name and directory if requested
  if ( Jali_BIN_OUTPUT_NAME )
    set_target_properties( ${target} PROPERTIES OUTPUT_NAME ${Jali_BIN_OUTPUT_NAME})
  endif()

  if ( Jali_BIN_OUTPUT_DIRECTORY )
    set_target_properties(${target} PROPERTIES OUTPUT_DIRECTORY ${Jali_BIN_OUTPUT_NAME})
  endif()



endfunction(ADD_Jali_EXECUTABLE)

#
# Jali Version Information:
# 
# Information about the current source is extracted from the mercurial repository and used to 
# create the version string (Jali_VERSION).  
#
# NOTE: this information won't be accessible without the full repository.
#       So for releases we need to extract this and set it as part of the tarball creation.
#
#   * if Jali_version.hh does not exist create it
#       * if mercurial is found
#            use mercurial to create version strings 
#       * else
#            use statically defined version strings
#       * endif
#   * endif
#   install Jali_version.hh
#

include(MercurialMacros)
include(PrintVariable)
include(InstallManager)

message(STATUS ">>>>>>>> JaliVersion.cmake")

if ( EXISTS ${CMAKE_SOURCE_DIR}/.hg/ )

find_package(Mercurial)
if ( NOT MERCURIAL_FOUND ) 

  message(ERROR "Could not locate Mercurial executable. Setting static version information")

  #
  # Not sure what is best for static information.  
  #
  set(Jali_VERSION_MAJOR 0)
  set(Jali_VERSION_MINOR 84)
  set(Jali_VERSION_PATCH dev)
  set(Jali_VERSION_HASH )

  #
  # Jali version
  #
  set(Jali_VERSION ${Jali_VERSION_MAJOR}.${Jali_VERSION_MINOR}-${Jali_VERSION_PATCH}_${Jali_VERSION_HASH})

else()

  mercurial_branch(Jali_HG_BRANCH)
  mercurial_global_id(Jali_HG_GLOBAL_HASH)
  mercurial_local_id(Jali_HG_LOCAL_ID)

  set(MERCURIAL_ARGS head ${Jali_HG_BRANCH} --template {latesttag}\n)
  execute_process(COMMAND  ${MERCURIAL_EXECUTABLE} ${MERCURIAL_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE Jali_HG_LATEST_TAG
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)


   # If Jali_HG_LATEST_TAG has a newline in it, chop it at the newline.
   STRING(FIND ${Jali_HG_LATEST_TAG} "\n" NEWLINE_INDEX)
   STRING(SUBSTRING ${Jali_HG_LATEST_TAG} 0 ${NEWLINE_INDEX} Jali_HG_LATEST_TAG)

   STRING(REGEX REPLACE "Jali-" "" Jali_HG_LATEST_TAG_VER ${Jali_HG_LATEST_TAG})	
   STRING(REGEX REPLACE "\\..*" "" Jali_HG_LATEST_TAG_MAJOR ${Jali_HG_LATEST_TAG_VER})	
   STRING(REGEX MATCH "\\.[0-9][0-9][\\.,-]" Jali_HG_LATEST_TAG_MINOR ${Jali_HG_LATEST_TAG_VER})	
   STRING(REGEX REPLACE "[\\.,-]" "" Jali_HG_LATEST_TAG_MINOR ${Jali_HG_LATEST_TAG_MINOR} )	

   set(Jali_VERSION_MAJOR ${Jali_HG_LATEST_TAG_MAJOR})
   set(Jali_VERSION_MINOR ${Jali_HG_LATEST_TAG_MINOR})

   #
   # Jali version
   #
   set(Jali_VERSION ${Jali_HG_LATEST_TAG_VER}_${Jali_HG_GLOBAL_HASH})

   STRING(REGEX REPLACE ".*\\.[0-9][0-9][\\.,-]" "" Jali_VERSION_PATCH ${Jali_VERSION})
   STRING(REGEX REPLACE ".*_" "" Jali_VERSION_HASH ${Jali_VERSION_PATCH})
   STRING(REGEX REPLACE "_.*" "" Jali_VERSION_PATCH ${Jali_VERSION_PATCH})

else()

  #
  # Not sure what is best for static information.  
  #
  set(Jali_VERSION_MAJOR 0)
  set(Jali_VERSION_MINOR 84)
  set(Jali_VERSION_PATCH dev)
  set(Jali_VERSION_HASH )

  #
  # Jali version
  #
  set(Jali_VERSION ${Jali_VERSION_MAJOR}.${Jali_VERSION_MINOR}-${Jali_VERSION_PATCH}_${Jali_VERSION_HASH})

endif()

# Status output
#message(STATUS "Jali Version:\t${Jali_VERSION}")
#message(STATUS "\tMercurial Branch:\t${Jali_HG_BRANCH}")
#message(STATUS "\tMercurial Global ID:\t${Jali_HG_GLOBAL_HASH}")
#message(STATUS "\tMercurial Local ID:\t${Jali_HG_LOCAL_ID}")
#message(STATUS "\tMercurial Tag Version:\t${Jali_HG_LATEST_TAG_VER}")


# Write the version header filed
set(version_template ${Jali_SOURCE_TOOLS_DIR}/cmake/Jali_version.hh.in)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/Jali_version.hh
               @ONLY)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/extras/Jali_version.hh
               @ONLY)

add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/Jali_version.hh)             

else()

 configure_file(${CMAKE_SOURCE_DIR}/Jali_version.hh
                ${CMAKE_CURRENT_BINARY_DIR}/Jali_version.hh
                @ONLY)
 add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/Jali_version.hh) 

 # extract the version from the include file
 file (STRINGS "${CMAKE_SOURCE_DIR}/Jali_version.hh" Jali_VERSION_HH REGEX "Jali_VERSION ")
 STRING(REGEX REPLACE ".*Jali_VERSION " "" Jali_VERSION ${Jali_VERSION_HH})
 STRING(STRIP ${Jali_VERSION} Jali_VERSION)

 STRING(REGEX REPLACE "\\..*" "" Jali_VERSION_MAJOR ${Jali_VERSION})
 STRING(REGEX MATCH "\\.[0-9][0-9]\\." Jali_VERSION_MINOR ${Jali_VERSION})
 STRING(REGEX REPLACE "\\." "" Jali_VERSION_MINOR ${Jali_VERSION_MINOR})
 STRING(REGEX REPLACE ".*\\..*\\." "" Jali_VERSION_PATCH ${Jali_VERSION})
 STRING(REGEX REPLACE ".*_" "" Jali_VERSION_HASH ${Jali_VERSION_PATCH})
 STRING(REGEX REPLACE "_.*" "" Jali_VERSION_PATCH ${Jali_VERSION_PATCH})

endif()

message(STATUS "\t >>>>>  Jali Version: ${Jali_VERSION}")
message(STATUS "\t >>>>>  MAJOR ${Jali_VERSION_MAJOR}")
message(STATUS "\t >>>>>  MINOR ${Jali_VERSION_MINOR}")
message(STATUS "\t >>>>>  PATCH ${Jali_VERSION_PATCH}")
message(STATUS "\t >>>>>  HASH  ${Jali_VERSION_HASH}")

