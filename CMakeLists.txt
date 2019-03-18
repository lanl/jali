# Copyright (c) 2019, Triad National Security, LLC
# All rights reserved.
#
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
# Jali 
#       Root CMakeLists.txt file
#

cmake_minimum_required(VERSION 3.8.0)
set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND})  

#if ("${CMAKE_PATCH_VERSION}" GREATER 11)
#  cmake_policy(SET CMP0022 NEW)
#endif()

#if ("${CMAKE_MAJOR_VERSION}" GREATER 2)
#  cmake_policy(SET CMP0037 OLD)
#  cmake_policy(SET CMP0048 OLD)
#endif()

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


# Some defines to indicate which partitioners are enabled
if (Jali_TPL_METIS_ENABLED)
  message(STATUS "Adding definition Jali_HAVE_METIS")
  add_definitions("-DJali_HAVE_METIS")
endif(Jali_TPL_METIS_ENABLED)
if (Jali_TPL_Zoltan_ENABLED)
  message(STATUS "Adding definition Jali_HAVE_ZOLTAN")
  add_definitions("-DJali_HAVE_ZOLTAN")
endif(Jali_TPL_Zoltan_ENABLED)
if ((NOT Jali_TPL_METIS_ENABLED) AND (NOT Jali_TPL_Zoltan_ENABLED))
  message(STATUS "No partitioner enabled?")
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


# compile examples only if requested

add_subdirectory(examples)
