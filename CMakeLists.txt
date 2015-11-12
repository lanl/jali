# -*- mode: cmake -*-

#
# Jali 
#       Root CMakeLists.txt file
#

# Require cmake 2.8.8 or higher
cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND})  

if ("${CMAKE_PATCH_VERSION}" GREATER 11)
  cmake_policy(SET CMP0022 NEW)
endif()

if ("${CMAKE_MAJOR_VERSION}" GREATER 2)
  cmake_policy(SET CMP0037 OLD)
  cmake_policy(SET CMP0048 OLD)
endif()

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# Define the project name
# This command will define
#     Jali_SOURCE_DIR
#     Jali_BINARY_DIR
project(Jali)

# Useful variables pointing to directories in the source tree
set(Jali_SOURCE_SRC_DIR      "${Jali_SOURCE_DIR}/src")
set(Jali_SOURCE_TOOLS_DIR    "${Jali_SOURCE_DIR}/tools")
set(Jali_SOURCE_CONFIG_DIR    "${Jali_SOURCE_DIR}/config")
set(Jali_PYTHON_DIR          "${Jali_SOURCE_DIR}/tools/py_lib")

# Set the module search path so find_package and include commands
# can locate files in <root source tree>/tools/cmake
set(Jali_MODULE_PATH "${Jali_SOURCE_CONFIG_DIR}/cmake")
set(CMAKE_MODULE_PATH 
    ${Jali_MODULE_PATH}
    ${Jali_MODULE_PATH}/Modules
    ${Jali_MODULE_PATH}/Utils)

# Code version defined
#include(JaliVersion)

# Addtional build options
include(JaliOptions)

# Find required Jali TPL
include(JaliTPL)

# The following line allows us to link third-party libraries not explicitly 
# found by the logic in JaliTPL.cmake.
link_directories(${CMAKE_INSTALL_PREFIX}/tpls/lib)

# Check the mesh framework choice
if ( (ENABLE_Unstructured)  AND
     (NOT ENABLE_STK_Mesh)  AND
     (NOT ENABLE_MOAB_Mesh) AND
     (NOT ENABLE_MSTK_Mesh)
   )
    message(FATAL_ERROR "Missing a mesh framework\n"
                        "Please enable at least one of the following mesh frameworks\n"
                        "-D ENABLE_MOAB_Mesh:BOOL=ON\n"
                        "-D ENABLE_STK_Mesh:BOOL=ON\n"
                        "-D ENABLE_MSTK_Mesh:BOOL=ON\n")
endif()    


# A property for accumulating the a global amanzi link line.
set_property(GLOBAL PROPERTY Jali_LINK_LINE "-L${CMAKE_INSTALL_PREFIX}/lib")

# A property for accumulating amanzi library targets
set_property(GLOBAL PROPERTY Jali_LIBRARY_TARGETS)
set(Jali_LINK_LINE_FILE "${Jali_BINARY_DIR}/link_line")  # A filename to write link-line to.
include(InstallManager)

# Source files for all binaries and libraries found under src
add_subdirectory(src)

# Define a list of all enabled TPLs. Must do this AFTER all the 
# CMakelist.txt files have been processed!
get_property(Jali_ENABLED_TPLS GLOBAL PROPERTY PACKAGES_FOUND)

# Make the exports only _after_ doing the build
create_exports()

option(ENABLE_Config_Report "Print out detailed information at the end of a configuration")
set(Jali_CONFIG_LOG "${Jali_BINARY_DIR}/jali-config.log"
    CACHE string "Jali configuration log file")
include(JaliConfigReport)
#if (ENABLE_Config_Report)
#    set(cat_exec "cat")
#    if (WIN32)
#        if( NOT UNIX)
#            set(cat_exec "type")
#        endif(NOT UNIX)
#    endif(WIN32)
#
#    execute_process(COMMAND "${cat_exec}" "${Jali_CONFIG_LOG}" OUTPUT_VARIABLE config_output)
#    print_variable(cat_exec)
#    print_variable(Jali_CONFIG_LOG)
#    print_variable(config_output)
#    message(STATUS "********************************************************************************")
#    message(STATUS "begin configuration output --\n${config_output}")
#    message(STATUS "end configuration output --")
#    message(STATUS "********************************************************************************")
#
#endif() 



# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
