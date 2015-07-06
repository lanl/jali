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

