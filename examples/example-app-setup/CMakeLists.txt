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


project (JALIAPP)

cmake_minimum_required(VERSION 2.8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#
# Jali Libraries and includes
#

if (NOT Jali_DIR)
  message(FATAL_ERROR "Error: Jali top level installation dir must be defined")
endif()

# Look for the Jali package

find_package(Jali)
if (Jali_FOUND) 
   message(STATUS "Located Jali")
   message(STATUS "Jali_DIR=${Jali_DIR}")
   message(STATUS "Jali_INCLUDE_DIRS=${Jali_INCLUDE_DIRS}")
   message(STATUS "Jali_LIBRARY_DIRS=${Jali_LIBRARY_DIRS}")
   message(STATUS "Jali_LIBRARIES=${Jali_LIBRARIES}")
   message(STATUS "Jali_TPL_INCLUDE_DIRS=${Jali_TPL_INCLUDE_DIRS}")
   message(STATUS "Jali_TPL_LIBRARY_DIRS=${Jali_TPL_LIBRARY_DIRS}")
   message(STATUS "Jali_TPL_LIBRARIES=${Jali_TPL_LIBRARIES}")
endif()


# Default build type is Release (Optimized)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()


# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${PROJECT_SOURCE_DIR})
endif ()

# if its an older CMAKE we have to explicitly enable C++11 support
if ("${CMAKE_MAJOR_VERSION}" LESS 3 OR "{CMAKE_PATCH_VERSION}" LESS 0)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif()


# include files and source files  
list(REVERSE Jali_LIBRARIES)

foreach (Jali_library ${Jali_LIBRARIES})
  message(STATUS "Adding Jali_library ${Jali_library} as an imported library")
  add_library(${Jali_library} STATIC IMPORTED)
  set_property(TARGET ${Jali_library} PROPERTY IMPORTED_LOCATION ${Jali_LIBRARY_DIRS}/lib${Jali_library}.a)
endforeach()

add_executable(testmesh src/main.cc)
include_directories(${PROJECT_SOURCE_DIR}/include ${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
target_link_libraries(testmesh ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

if ("${CMAKE_MAJOR_VERSION}" GREATER 2 AND "{CMAKE_PATCH_VERSION}" GREATER 0)
  set_property(TARGET testmesh PROPERTY CXX_STANDARD 11)
  set_property(TARGET testmesh PROPERTY CXX_STANDARD_REQUIRED ON)
endif()

