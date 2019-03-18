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
# Functions for building a link-line for linking to Jali.
#
# See the modified CMakeLists.txt files for usage. 
#
# This functionality should be merged with InstallManager.

include(ParseLibraryList)
include(PrintVariable)

# From JaliConfigReport.cmake:
set(build_timestamp "Not available on this platform")
if (UNIX)
  execute_process(COMMAND "date"
    RESULT_VARIABLE _ret_code
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE  _stderr
    )
  string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()

function(_parse_add_libraries library_list libraries_to_add)

if (library_list)
  parse_library_list(
    ${library_list}
    FOUND   libraries_split
    DEBUG   debug_libraries
    OPT     opt_libraries
    GENERAL general_libraries)

  if (libraries_split)
    message("Libraries for ${package} are present in multiple debug and/or opt versions")
    if(${CMAKE_BUILD_TYPE} MATCHES "debug")
      message("Adding debug libraries")
      set(${libraries_to_add} "${debug_libraries}" PARENT_SCOPE)
    else()
      message("Adding optimized libraries")
      set(${libraries_to_add} "${opt_libraries}" PARENT_SCOPE)
    endif()
  else()
    set(${libraries_to_add} "${library_list}" PARENT_SCOPE)
  endif()

endif()
endfunction()





macro(_add_to_link_line)
  set_property(GLOBAL APPEND PROPERTY Jali_LINK_LINE ${ARGV})
endmacro()

macro(_add_to_target_list)
  set_property(GLOBAL APPEND PROPERTY Jali_LIBRARY_TARGETS ${ARGV})
endmacro()



macro(_add_directories_to_link_line)
  foreach(directory ${ARGV})
    _add_to_link_line("-L${directory}")
  endforeach()
endmacro()

macro(_add_libraries_to_link_line)
  foreach(library ${ARGV})
    if(EXISTS ${library})  
      _add_to_link_line("${library}")   # If it's a filename, add it as given.
    else()
      _add_to_link_line("-l${library}") # Else, add it as a library to be looked up.
    endif()
  endforeach()
endmacro()



macro(add_package_libraries)

  # Grab the project name to find the dependent libraries
  SET(package ${PROJECT_NAME})

  # Add the directory locations of libraries it depends on.
  _add_directories_to_link_line("${${package}_LIBRARY_DIR}")
  _add_directories_to_link_line("${${package}_LIBRARY_DIRS}")

  # ${package}_LIBRARIES may contain debug and opt keywords, so parse the list into to_add:
  _parse_add_libraries("${${package}_LIBRARIES}" to_add)

  _add_libraries_to_link_line("${to_add}")

  # Add the accumulated contents of value_list to the link line property
  set_property(GLOBAL APPEND PROPERTY Jali_LINK_LINE ${value_list})

endmacro()



macro(add_Jali_libraries libraries)
  _add_libraries_to_link_line(${libraries})
  _add_to_target_list(${libraries})   
endmacro()

