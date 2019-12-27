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
#   software without specific prior written permission.
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
# Jali Build Options
#
#
# This file is intended define build options
# related to compile options, build types, etc.
# Options related to Third Party Libraries (TPL)
# can be found in JaliTPL.cmake

# Standard CMake modules
include(CMakeDependentOption)
include(FeatureSummary)

enable_language(C)
enable_language(CXX)

# ENABLE C++11 support

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()


## No idea why we need this.
## I think it was required for Franklin build. -- lpritch
#if(PREFER_STATIC_LIBRARIES)
#  # Prefer static libraries, but don't require that everything must be static. 
#  set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
#endif(PREFER_STATIC_LIBRARIES)

#if(BUILD_STATIC_EXECUTABLES)
#    set(CMAKE_EXE_LINKER_FLAGS -static)
#    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
#    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
#    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
#    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
#    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
#    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
##    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
#endif(BUILD_STATIC_EXECUTABLES)

#
# Options
# 
# Testing
cmake_dependent_option(ENABLE_TESTS "Enable unit testing" ON
                       "ENABLE_UnitTest" ON)
add_feature_info(TESTS
                 ENABLE_TESTS
                 "Toggle for unit tests")
if (ENABLE_TESTS)
    set(BUILD_TESTS 1)
endif()    

# Some platforms require all binaries linking to MPI
# only run through the MPIEXEC binary
option(TESTS_REQUIRE_MPIEXEC "Run all tests with the MPIEXEC binary" FALSE)

# Need this option if the PATH environment does not include '.'
option(TESTS_REQUIRE_FULLPATH "Append full path to test binaries" TRUE)

