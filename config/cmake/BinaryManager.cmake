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

