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
#  Jali
#    Error Handling, Throw, Catch, etc.
#

# Jali module, include files found in JALI_MODULE_PATH
include(TestManager)
#include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   ERROR_SOURCE_DIR
#   ERROR_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ERROR)

set(ERROR_headers
  errors.hh
  exceptions.hh)
list(TRANSFORM ERROR_headers PREPEND "${ERROR_SOURCE_DIR}/")

set(ERROR_sources
  errors.cc
  exceptions.cc)


#
# Library Target: error_handling
#
add_library(jali_error_handling ${ERROR_sources})
set_target_properties(jali_error_handling PROPERTIES PUBLIC_HEADER "${ERROR_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_error_handling ALIAS jali_error_handling)

# We should change install_interface to include/error_handling but that requires
# changing all client codes
target_include_directories(jali_error_handling PUBLIC
  $<BUILD_INTERFACE:${ERROR_BINARY_DIR}>
  $<BUILD_INTERFACE:${ERROR_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (BUILD_TESTS)

  # Create the error test suite
  file(GLOB error_test_files
    test/test_errors.cc
    test/test_exceptions.cc
    )

  message(STATUS "UnitTest++_LIBRARIES ${UnitTest++_LIBRARIES}")
  add_Jali_test(test_error_handling test_error_handling
    KIND unit
    SOURCE test/Main.cc ${error_test_files}
    LINK_LIBS jali_error_handling ${UnitTest++_LIBRARIES})

endif()



