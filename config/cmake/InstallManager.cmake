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
# Functions for managing the install targets
#


include(CMakeParseArguments)
include(JaliLinkLine)
include(CreateImportedTargetFile)


export(PACKAGE Jali)

#
# Usage: ADD_INSTALL_INCLUDE_FILE( file1 file2 file3 ... )
#
# Arguments:
#  A list of files that will be installed in the Jali_INSTALL_INCLUDE_DIR
#
#
function ( ADD_INSTALL_INCLUDE_FILE )

  foreach(_inc_file ${ARGV})
    install(
      FILES ${_inc_file}
      DESTINATION include
      )
  endforeach()

endfunction( ADD_INSTALL_INCLUDE_FILE )

#
# Usage: ADD_INSTALL_LIBRARY( lib1 lib2 lib3 ... )
#
# Arguments:
#  A list of libraries that will be installed in the Jali_INSTALL_LIB_DIR
#
#
function ( ADD_INSTALL_LIBRARY )

  install(
    TARGETS ${ARGV}
    EXPORT JaliTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  # Add the libraries to our global list
  add_Jali_libraries(${ARGV})

  # Add dependency libaries as determined by the package definition.
  add_package_libraries()

endfunction( ADD_INSTALL_LIBRARY )


#
# Usage: ADD_INSTALL_SHELLSCRIPT( script1 script2 ... )
#
# Arguments:
#  A list of shell scripts that will be installed in the Jali_INSTALL_BIN_DIR
#
#
function ( ADD_INSTALL_SHELLSCRIPT )

foreach(_shellscript_file ${ARGV})
  install(
    FILES ${_shellscript_file}
    DESTINATION bin
    )
endforeach()

endfunction( ADD_INSTALL_SHELLSCRIPT )




#
# Usage: ADD_INSTALL_BINARY( exe1 exe2 ... )
#
# Arguments:
#  A list of executables that will be installed in the Jali_INSTALL_BIN_DIR
#
#
function ( ADD_INSTALL_BINARY )

foreach(_bin_file ${ARGV})
  install(
    TARGETS ${_bin_file}
    EXPORT JaliTargets
    DESTINATION bin
    )
endforeach()

endfunction( ADD_INSTALL_BINARY )

#
# Usage: create_tpl_export_file( <package list> | package1 package2 ... )
#
# Arguments: Semicolon deliminated list of TPL package names or
#            individual package names
#
function( CREATE_TPL_EXPORT_FILE )

  # Print the usage
  macro( _print_usage )
    message("\nUsage: create_tpl_export_file( outfile PACKAGES <package list> | package1 package2 ... )\n")
  endmacro()

  # Parse Arguments
  set(_options "")
  set(_oneValue "")
  set(_multiValue "PACKAGES")
  cmake_parse_arguments(BUILD_TPL "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN})

  if (NOT BUILD_TPL_PACKAGES)
    _print_usage()
    message(FATAL_ERROR "Require a package list to build export file")
  endif()

  list(GET BUILD_TPL_UNPARSED_ARGUMENTS 0 BUILD_TPL_OUTFILE)
  if (NOT BUILD_TPL_OUTFILE)
    _print_usage()
    message(FATAL_ERROR "Must define an output file")
  endif()  

  # BEGIN MACROS

  # Write the header for the file
  macro(_write_header)

    file(WRITE ${BUILD_TPL_OUTFILE} "# ############################################################################ #\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# Jali TPL (External Software Packages) Configuration\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# ############################################################################ #\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "\n")

    file(APPEND ${BUILD_TPL_OUTFILE} "# ------------------------------------------------------------------------------\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# TPL Config File Directory\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# ------------------------------------------------------------------------------\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "GET_FILENAME_COMPONENT(SELF_DIR \${CMAKE_CURRENT_LIST_FILE} PATH)\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "\n")

    file(APPEND ${BUILD_TPL_OUTFILE} "# ------------------------------------------------------------------------------\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# Imported Target File\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# ------------------------------------------------------------------------------\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "SET(Jali_IMPORT_TARGET_FILE \${SELF_DIR}/JaliImportedTargets.cmake)\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "INCLUDE(\${Jali_IMPORT_TARGET_FILE})\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "\n")

  endmacro(_write_header)

  # Write a CMake set variable command
  macro(_write_cmake_variable _cmake_var_name _var_value)

    if (${_var_value})
      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_cmake_var_name} ${${_var_value}})\n")
    endif()  

  endmacro(_write_cmake_variable)  

  # Add package to the file
  macro( _add_package _package )

    file(APPEND ${BUILD_TPL_OUTFILE} "\n#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# TPL: ${_package}\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")

    set(found_package_flag ${_package}_FOUND)
    set(_package_enabled_flag_var "Jali_TPL_${_package}_ENABLED")
    set(_package_dir_var          "Jali_TPL_${_package}_DIR")

    if ( ${found_package_flag} ) 

      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_enabled_flag_var} ON)\n")
      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_dir_var} ${${_package}_DIR})\n")

      set(_var_name_list "INCLUDE_DIR;INCLUDE_DIRS;LIBRARIES;LIBRARY_DIR;LIBRARY_DIRS")
      foreach (_append_name ${_var_name_list})
	set(_write_var_name "Jali_TPL_${_package}_${_append_name}")
	set(_var_name       "${_package}_${_append_name}")
	_write_cmake_variable(${_write_var_name} ${_var_name})
      endforeach(_append_name)	

#
#     Also include CMake configuration files for the packages so that 
#     any recursive dependencies are discovered for linking
#

      find_file(${_package}_config_file NAMES ${_package}Config.cmake PATHS ${${_package}_DIR} NO_DEFAULT_PATH)
      if (${_package}_config_file-NOTFOUND)
	find_file(${_package}_config_file NAMES ${_package}-config.cmake PATHS ${${_package}_DIR} NO_DEFAULT_PATH)
      endif()

      if (${_package}_config_file)
#	message(status "${_package} config dir ------ ${${_package}_DIR}")
#	message(status "${_package} config file ------ ${${_package}_config_file}")
	file(APPEND ${BUILD_TPL_OUTFILE} "include(\"${${_package}_config_file}\")")
      endif()

    else()

      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_enabled_flag_var} OFF)\n")

    endif()

  endmacro(_add_package) 


  # Begin MAIN

  if (NOT EXISTS ${BUILD_TPL_OUTFILE})
    _write_header(${BUILD_TPL_OUTFILE})
  endif()  

  # Loop through each package and update the output file
  foreach(_package IN LISTS BUILD_TPL_PACKAGES)
    #print_variable(_package)
    _add_package(${_package}) 
  endforeach() 

  # Also, at the end build a concatenated list of all the enabled TPL
  # libaries, their include dirs and their library files. Should probably
  # do this inside _add_package but I don't know how to

  foreach(_package IN LISTS BUILD_TPL_PACKAGES)
    if (${${_package}_FOUND})
      list(APPEND Jali_TPL_LIBRARIES ${${_package}_LIBRARIES})
      list(APPEND Jali_TPL_LIBRARY_DIRS ${${_package}_LIBRARY_DIR})
      list(APPEND Jali_TPL_LIBRARY_DIRS ${${_package}_LIBRARY_DIRS})
      list(APPEND Jali_TPL_INCLUDE_DIRS ${${_package}_INCLUDE_DIR})
      list(APPEND Jali_TPL_INCLUDE_DIRS ${${_package}_INCLUDE_DIRS})
      set(tpl_libraries    "${${package}_LIBRARIES}")
      foreach( extra_tpl_library ${tpl_libraries} )
	get_filename_component(extra_library_path ${extra_tpl_library} PATH)
	list(APPEND Jali_TPL_LIBRARY_DIRS ${extra_library_path})
      endforeach()
    endif ()
  endforeach()
  list(REMOVE_DUPLICATES Jali_TPL_LIBRARIES)
  list(REMOVE_DUPLICATES Jali_TPL_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES Jali_TPL_INCLUDE_DIRS)

  # write the concatenated list variables out

  file(APPEND ${BUILD_TPL_OUTFILE} "\n#\n# All TPL LIBRARIES\n#\n")
  file(APPEND ${BUILD_TPL_OUTFILE} "set(Jali_TPL_LIBRARIES ${Jali_TPL_LIBRARIES})\n")

  file(APPEND ${BUILD_TPL_OUTFILE} "\n#\n# All TPL LIBRARY DIRECTORIES\n#\n")
  file(APPEND ${BUILD_TPL_OUTFILE} "set(Jali_TPL_LIBRARY_DIRS ${Jali_TPL_LIBRARY_DIRS})\n")

  file(APPEND ${BUILD_TPL_OUTFILE} "\n#\n# All TPL INCLUDE DIRECTORIES\n#\n")
  file(APPEND ${BUILD_TPL_OUTFILE} "set(Jali_TPL_INCLUDE_DIRS ${Jali_TPL_INCLUDE_DIRS})\n")
  message(STATUS "Leaving create_tpl_export_file") 

endfunction ( CREATE_TPL_EXPORT_FILE )

#
# Usage: makefile_include_dirs(CMAKE_INCLUDE_LIST in_list 
#                               MAKE_INCLUDE_LIST out_list)
# 
# Arguments:
#          CMAKE_INCLUDE_LIST List of include directories
#          MAKE_INCLUDE_LIST  List include directories for make
# 
function(makefile_include_dirs)

    cmake_parse_arguments(PARSE_ARGS "" "MAKE_INCLUDE_LIST" "CMAKE_INCLUDE_LIST" ${ARGN})
    #print_variable(PARSE_ARGS_CMAKE_INCLUDE_LIST)
    #print_variable(PARSE_ARGS_MAKE_INCLUDE_LIST)

    set(tmp_inc_list)
    set(loop_list ${PARSE_ARGS_CMAKE_INCLUDE_LIST})
    list(REMOVE_DUPLICATES loop_list)
    foreach( dir  ${loop_list})
      set(i_path "-I${dir} ")
      list(APPEND tmp_inc_list ${i_path})
    endforeach() 

    set(tmp_make_list)
    string(REGEX REPLACE ";" "" tmp_make_list ${tmp_inc_list})
    set(${PARSE_ARGS_MAKE_INCLUDE_LIST} "${tmp_make_list}" PARENT_SCOPE)

endfunction(makefile_include_dirs)

#
# Usage: makefile_library_dirs(CMAKE_LIB_LIST in_list 
#                              MAKE_LIB_LIST out_list)
# 
# Arguments:
#          CMAKE_LIB_LIST List of library directories
#          MAKE_LIB_LIST  List library directories for make
# 
function(makefile_library_dirs)

    cmake_parse_arguments(PARSE_ARGS "" "MAKE_LIB_LIST" "CMAKE_LIB_LIST" ${ARGN})
    #print_variable(PARSE_ARGS_CMAKE_LIB_LIST)
    #print_variable(PARSE_ARGS_MAKE_LIB_LIST)

    set(tmp_lib_list)
    set(loop_list ${PARSE_ARGS_CMAKE_LIB_LIST})
    list(REVERSE loop_list)
    list(REMOVE_DUPLICATES loop_list)
    list(REVERSE loop_list)
    foreach( dir  ${loop_list})
      set(l_path "-L${dir} ")
      list(APPEND tmp_lib_list ${l_path})
    endforeach() 

    set(tmp_make_list)
    string(REGEX REPLACE ";" "" tmp_make_list ${tmp_lib_list})
    set(${PARSE_ARGS_MAKE_LIB_LIST} "${tmp_make_list}" PARENT_SCOPE)

endfunction(makefile_library_dirs)

#
# Usage: create_exports
#
# Arguments: None
#
#
function (CREATE_EXPORTS)

# Template file located in the CMake module directory

# Define Jali_INCLUDE_DIRS, Jali_LIBRARY_DIRS, Jali_LIBRARIES
set(Jali_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
set(Jali_LIBRARY_DIRS "${CMAKE_INSTALL_PREFIX}/lib")
get_property(Jali_LIBRARIES GLOBAL PROPERTY Jali_LIBRARY_TARGETS)
list(REVERSE Jali_LIBRARIES)

# Find the packages found for Jali
get_property(Jali_TPL_LIST GLOBAL PROPERTY PACKAGES_FOUND)
get_property(LINK_LINE GLOBAL PROPERTY Jali_LINK_LINE)

# Define Jali_TPL_INCLUDE_DIRS, Jali_TPL_LIBRARY_DIRS, Jali_TPL_LIBRARIES

#foreach( package ${Jali_TPL_LIST} )
#  set(tpl_include_dir "${${package}_INCLUDE_DIR}")
#  set(tpl_include_dirs "${${package}_INCLUDE_DIRS}")
#  list(APPEND Jali_TPL_INCLUDE_DIRS ${tpl_include_dir} ${tpl_include_dirs})
#
#  set(tpl_library_dir  "${${package}_LIBRARY_DIR}")
#  set(tpl_library_dirs "${${package}_LIBRARY_DIRS}")
#  list(APPEND Jali_TPL_LIBRARY_DIRS ${tpl_library_dir} ${tpl_library_dirs})
#
#  set(tpl_libraries    "${${package}_LIBRARIES}")
#  foreach( extra_tpl_library ${tpl_libraries} )
#    get_filename_component(extra_library_path ${extra_tpl_library} PATH)
#    list(APPEND Jali_TPL_LIBRARY_DIR ${extra_library_path})
#  endforeach()
#endforeach()
#list(REMOVE_DUPLICATES Jali_TPL_INCLUDE_DIRS)
#list(REMOVE_DUPLICATES Jali_TPL_LIBRARY_DIRS)

# Convert the link line to a space deliminated string
foreach (arg ${LINK_LINE})
  set(LINK_LINE_STRING "${LINK_LINE_STRING} ${arg}")
endforeach()

# Write and install the link-line file
file(WRITE ${Jali_LINK_LINE_FILE} ${LINK_LINE_STRING})
install(FILES ${Jali_LINK_LINE_FILE} DESTINATION lib)

# Write the imported target file
set(import_target_file "${Jali_BINARY_DIR}/JaliImportedTargets.cmake")
create_imported_target_file(${import_target_file})
install(FILES ${import_target_file} DESTINATION lib)

# Write the TPL file
set(tpl_config_file "${Jali_BINARY_DIR}/JaliConfigTPL.cmake")
if ( EXISTS ${tpl_config_file} )
  file(REMOVE ${tpl_config_file})
endif()  
create_tpl_export_file(${tpl_config_file}
                       PACKAGES ${Jali_ENABLED_TPLS})
install(FILES ${tpl_config_file} DESTINATION lib)				   

# Write the export Makefile and add to the include install list
makefile_include_dirs(CMAKE_INCLUDE_LIST ${Jali_INCLUDE_DIRS}
                      MAKE_INCLUDE_LIST Jali_MAKE_INCLUDE_DIRS) 
makefile_library_dirs(CMAKE_LIB_LIST ${Jali_LIBRARY_DIRS}
                      MAKE_LIB_LIST Jali_MAKE_LIBRARY_DIRS) 
set(in_makefile  "${Jali_MODULE_PATH}/MakefileConfig.export.in")
set(out_makefile "${Jali_BINARY_DIR}/Makefile.export")
configure_file("${in_makefile}" "${out_makefile}")
install(FILES "${out_makefile}" DESTINATION lib)

# Write the JaliConfig.cmake file
set(in_config   "${Jali_MODULE_PATH}/JaliConfig-install.cmake.in")
set(out_config   "${Jali_BINARY_DIR}/JaliConfig.cmake")
configure_file(${in_config} ${out_config})
install(FILES ${out_config} DESTINATION lib)

# Write the JaliConfigVersion.cmake file
#set(in_config   "${Jali_MODULE_PATH}/JaliConfigVersion-install.cmake.in")
#set(out_config   "${Jali_BINARY_DIR}/JaliConfigVersion.cmake")
#configure_file(${in_config} ${out_config} @ONLY)
#install(FILES ${out_config} DESTINATION lib)

# Write the CMake configuration target file
message(STATUS "Writing target file")
install(EXPORT JaliTargets
        DESTINATION lib
	NAMESPACE Jali_
	FILE JaliTargets.cmake)

# XML schema 
#install(FILES ${Jali_SOURCE_DIR}/doc/input_spec/schema/Jali.xsd DESTINATION bin)

# If MSTK utilities were found (this processing is in FindMSTK.cmake), install them in bin
if (MSTK_UTILITIES)
   install(PROGRAMS ${MSTK_UTILITIES} DESTINATION bin)
endif()

# Write the evaluator registration macro file
#message(STATUS "Writing evaluator registration macro file")
#install(FILES tools/cmake/RegisterEvaluators.cmake
#        DESTINATION lib)

endfunction()

