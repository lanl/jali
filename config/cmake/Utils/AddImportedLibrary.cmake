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
# Jali Print variable
#
#
# ADD_IMPORTED_LIBRARY(name [SHARED | STATIC]
#                      LOCATION <path>
#                      [ LINK_LANGUAGES <lang1> <lang2> <lang3> ... ]
#                      [ LINK_INTERFACE_LIBRARIES <lib1> <lib2> ... ]
#                     )
#                    
#                      

include(CMakeParseArguments)
include(PrintVariable)
function(ADD_IMPORTED_LIBRARY target_name)

  set(_options SHARED STATIC)
  set(_oneValueArgs LOCATION)
  set(_multiValueArgs LINK_LANGUAGES LINK_INTERFACE_LIBRARIES)

  cmake_parse_arguments(PARSE "${_options}" "${_oneValueArgs}" "${_multiValueArgs}" ${ARGN} ) 

  # --- Check what has been passed in

  # SHARED and STATIC can not be set at the same time
  if ( "${PARSE_STATIC}" AND "${PARSE_SHARED}" )
    message(FATAL_ERROR "Can not specify imported library as shared and static.")
  endif()

  # Require a location
  if ( NOT PARSE_LOCATION )
    message(FATAL_ERROR "Must specify a location to define an imported library target.")
  endif()

  # Check to see if name already exists as a target
  if ( TARGET "${target_name}" )
    message(FATAL_ERROR "Target ${name} is already defined.")
  endif()


  # --- Set the library type
  set(lib_type UNKNOWN)
  if(PARSE_STATIC)
    set(lib_type STATIC)
  endif()  
  if(PARSE_SHARED)
    set(lib_type SHARED)
  endif()

  # --- Add the library target 
  add_library(${target_name} ${lib_type} IMPORTED)

  # --- Update the global property that tracks imported targets
  set(prop_name IMPORTED_${lib_type}_LIBRARIES)
  get_property(prop_value GLOBAL PROPERTY ${prop_name})
  set_property(GLOBAL PROPERTY ${prop_name} ${prop_value} ${target_name})

  # --- Set the properties
  set_target_properties(${target_name} PROPERTIES
                        IMPORTED_LOCATION ${PARSE_LOCATION})
  if ( PARSE_LINK_LANGUAGES )
    set_target_properties(${target_name} PROPERTIES
                        IMPORTED_LINK_INTERFACE_LANGUAGES "${PARSE_LINK_LANGUAGES}")
  endif()
  if ( PARSE_LINK_INTERFACE_LIBRARIES )
    set_target_properties(${target_name} PROPERTIES
                          IMPORTED_LINK_INTERFACE_LIBRARIES "${PARSE_LINK_INTERFACE_LIBRARIES}")
  endif()
     
  
endfunction(ADD_IMPORTED_LIBRARY)

