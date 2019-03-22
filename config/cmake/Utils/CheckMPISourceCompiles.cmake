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
# CHECK_MPI_SOURCE_COMPILES
#
#  A function that verifies that the CMAKE_C_COMPILER, 
#  CMAKE_CXX_COMPILER and CMAKE_Fortran_COMPILER can compile and link 
#  a MPI test file. Since this function is a simple wrapper
#  call to check_c{xx}_source_compiles macros, the following
#  variables can be set before the this function to alter
#  how the mocros behave:
#         CMAKE_REQUIRED_FLAGS = string of compile command line flags
#         CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#         CMAKE_REQUIRED_INCLUDES = list of include directories
#         CMAKE_REQUIRED_LIBRARIES = list of libraries to link
#
function(CHECK_MPI_SOURCE_COMPILES FLAG)

  include(CheckCSourceCompiles)

  message(STATUS "Checking whether C compiler is MPI wrapper")

  set(_mpi_c_source "
      #include <stdio.h>
      #include <mpi.h>
      void main(int argc, char **argv) 
      {
        MPI_Init(&argc,&argv);
        puts(__LINE__);
        MPI_Finalize();
      }")

    # CHECK_C_SOURCE_COMPILES is a MACRO not a function
    # check $VAR for compile status.
    check_c_source_compiles("${_mpi_c_source}" MPI_C)
    if(${MPI_C})
      message(STATUS "Checking whether C compiler is MPI wrapper - yes")
    else()
      message(STATUS "Checking whether C compiler is MPI wrapper - no")
      set(MPI_C False CACHE INTERNAL "Test MPI_C")
    endif()   
        

    # MPI C++ bindings have disappeared in recent MPI versions - skip check
    include(CheckCXXSourceCompiles)
    set(_mpi_cxx_source "
        #include <iostream>
        #include \"mpi.h\"
        int main(int argc, char *argv[])
        {
          MPI_Init(argc,argv);
          MPI_Finalize();
          return 0;
        }")
    
    # CHECK_C_SOURCE_COMPILES is a MACRO not a function
    # check $VAR for compile status.
    check_cxx_source_compiles("${_mpi_cxx_source}" MPI_CXX)

    if(${MPI_CXX})
      message(STATUS "Checking whether CXX compiler is MPI wrapper - yes")
    else()
      message(STATUS "Checking whether CXX compiler is MPI wrapper - no")
      set(MPI_CXX False CACHE INTERNAL "Test MPI_CXX")
    endif() 

     # Now push up the flag to the parent namespace
     if ( ${MPI_C} AND ${MPI_CXX} )
       set(${FLAG} TRUE PARENT_SCOPE)
     else()
       set(${FLAG} FALSE PARENT_SCOPE)
     endif()
         
endfunction(CHECK_MPI_SOURCE_COMPILES)
