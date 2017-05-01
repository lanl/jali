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

  message(STATUS "Checking whether C compiler can compile MPI program")

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
      message(STATUS "Checking whether C compiler can compile MPI program - yes")
    else()
      message(STATUS "Checking whether C compiler can compile MPI program - no")
    endif()   
        


    include(CheckCXXSourceCompiles)
    message(STATUS "Checking whether C++ compiler can compile MPI program")
    set(_mpi_cxx_source "
        #include <iostream>
        #include \"mpi.h\"
        int main(int argc, char *argv[])
        {
          MPI::Init(argc,argv);
          MPI::Finalize();
          return 0;
        }")
    
    # CHECK_C_SOURCE_COMPILES is a MACRO not a function
    # check $VAR for compile status.
    check_cxx_source_compiles("${_mpi_cxx_source}" MPI_CXX)

    if(${MPI_CXX})
      message(STATUS "Checking whether C++ compiler can compile MPI program - yes")
    else()
      message(STATUS "Checking whether C++ compiler can compile MPI program - no")
    endif() 

    # Check the Fortran compiler if enabled
    if ( CMAKE_Fortran_COMPILER ) 
      message(STATUS "Checking whether Fortran compiler can compile MPI program")
      set(_mpi_fortran_source "
          PROGRAM mpi_test
          INTEGER ierr
          call MPI_Init(ierr)
          call MPI_Finalize(ierr)
          STOP
          END PROGRAM")

      # As of CMAKE 2.8.6 no CheckFortranSourceCompiles macro
      set(test_file "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_mpi.f") 
      file(WRITE "${test_file}" "${_mpi_fortran_source}\n")

      message(STATUS "Performing Test MPI_Fortran")
      try_compile(MPI_Fortran
                  ${CMAKE_BINARY_DIR}
                  ${test_file}
                  OUTPUT_VARIABLE result)
      if(${MPI_Fortran})
        message(STATUS "Performing Test MPI_Fortran - Success")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                    "Performing Fortran SOURCE FILE Test MPI_Fortran succeeded with the following output:\n"
                    "${result}\n"
                    "Source file was:\n${_mpi_fortran_source}\n")
        message(STATUS "Checking whether Fortran compiler can compile MPI program - yes")
      else()
        message(STATUS "Performing Test MPI_Fortran - Failed")
        set(MPI_Fortran "" CACHE INTERNAL "Test MPI_Fortran")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                    "Performing Fortran SOURCE FILE Test MPI_Fortran failed with the following output:\n"
                    "${result}\n"
                    "Source file was:\n${_mpi_fortran_source}\n")
        message(STATUS "Checking whether Fortran compiler can compile MPI program - no")
       endif()            

     endif()

     # Now push up the flag to the parent namespace
     if ( CMAKE_Fortran_COMPILER )
       if ( "${MPI_C}" AND "${MPI_CXX}" AND "${MPI_Fortran}" )
         set(${FLAG} TRUE PARENT_SCOPE)
       else()
         set(${FLAG} FALSE PARENT_SCOPE)
       endif()
     else()
       if ( "${MPI_C}" AND "${MPI_CXX}" )
         set(${FLAG} TRUE PARENT_SCOPE)
       else()
         set(${FLAG} FALSE PARENT_SCOPE)
       endif()
     endif()  
         
endfunction(CHECK_MPI_SOURCE_COMPILES)
