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
# Functions for managing tests.
#

include(CMakeParseArguments)
include(PrintVariable)

function(_APPEND_TEST_LABEL test_name label)

  get_test_property(${test_name} LABELS current_labels)
  if (current_labels)
    set_tests_properties(${test_name} PROPERTIES LABELS "${current_labels};${label}")
  else()  
    set_tests_properties(${test_name} PROPERTIES LABELS "${label}")
  endif()

endfunction(_APPEND_TEST_LABEL)

function(_ADD_TEST_KIND_LABEL test_name kind_in)

  set(kind_prefixes UNIT INT REG Jali)

  string(TOUPPER "${kind_in}" kind)

  foreach(kind_prefix ${kind_prefixes})
    string(REGEX MATCH "${kind_prefix}" match ${kind})
    if(match)
      break()
    endif()
  endforeach()

 if (match)
    _append_test_label(${test_name} ${match})
  else()
    message(FATAL_ERROR "Invalid test label ${kind_in} (Valid Labels:${kind_prefixes})")
  endif()

endfunction(_ADD_TEST_KIND_LABEL)


# Usage:
#
# ADD_Jali_TEST(<test_name> <test_executable>
#                  [arg1 ...]
#                  KIND [unit | int | reg | Jali ]
#                  [Jali_INPUT file.xml]
#                  [SOURCE file1 file2  ...]
#                  [LINK_LIBS lib1 lib2 ...]
#                  [DEPENDS test1 test2 ...]
#                  [PARALLEL] [EXPECTED_FAIL]
#                  [NPROCS procs1 ... ]
#                  [MPI_EXEC_ARGS arg1 ... ])

#
# Arguments:
#  test_name: the name given to the resulting test in test reports
#  test_executable: The test executable which performs the test. 
#                   Required if KIND is unit, int or reg
#  arg1 ...: Additional arguments for the test executable
#
# Keyword KIND is required and should be one of unit, int, reg or Jali.
#         Jali is a special case where the test executable is
#         set to the main Jali binary.
#
# KEYWORD Jali_INPUT is required if keyword KIND is set to Jali. This
#         key word defines the Jali XML input file.
#
# Option PARALLEL signifies that this is a parallel job. This is also
# implied by an NPROCS value > 1
#
# Optional NPROCS keyword starts a list of the number of processors to
# run the test on. Defaults to 1.
#
# Optional MPI_EXEC_ARGS keyword denotes extra arguments to give to
# mpi. It is ignored for serial tests.
#
# Optional SOURCE_FILES keyword that defines a list of source files
# required to build test_executable. An add_executable call will be made
# if this option is active.
#
# Optional LINK_LIBS keyword defines a list of link libraries or link options
# to link test_executable. An target_link_libraries call will be made if
# this option is active.
#
# Optional DEPENDS keyword defines a list of tests that should finish before
# test_name

function(ADD_Jali_TEST test_name)

  # --- Initialize 

  # Check test_name
  if ( NOT test_name )
    message(FATAL_ERROR "Must define a test name.")
  endif()

  # Parse through the remaining options
  set(options PARALLEL EXPECTED_FAIL)
  set(oneValueArgs KIND Jali_INPUT)
  set(multiValueArgs NPROCS SOURCE LINK_LIBS DEPENDS MPI_EXEC_ARGS)
  cmake_parse_arguments(Jali_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # --- Check options

  # Require a KIND value
  if ( NOT Jali_TEST_KIND )
    message(FATAL_ERROR "A test type has not been specified for ${test_name}.")
  endif()

  # Force each test to parallel run if mpiexec is required
  if(TESTS_REQUIRE_MPIEXEC)
    set(Jali_TEST_PARALLEL TRUE)
  endif()

  # Force each PARALLEL TRUE if NPROCS set 
  if(Jali_TEST_NPROCS AND ( "${Jali_TEST_NPROCS}" GREATER 1 ) )
    set(Jali_TEST_PARALLEL TRUE)
  endif()  

  # Default to nprocs=1 when running parallel
  if ( Jali_TEST_PARALLEL AND (NOT Jali_TEST_NPROCS) )
    set(Jali_TEST_NPROCS 1)
  endif() 

  # Test the value of number of procs value
  if(Jali_TEST_NPROCS)
    if(NOT ("${Jali_TEST_NPROCS}" GREATER 0) )
      message(FATAL_ERROR "${Jali_TEST_NPROCS} is an invalid NPROCS value.")
    endif()

    if(MPI_EXEC_MAX_NUMPROCS AND Jali_TEST_PARALLEL)
      if ( "${MPI_EXEC_MAX_NUMPROCS}" LESS "${Jali_TEST_NPROCS}")
        message(WARNING "Test ${test_name} request too many nprocs (${Jali_TEST_NPROCS}). "
                        "Will skip this test.")
        return()
      endif()
    endif()
  endif()

  # --- Define the test executable

  if ( "${Jali_TEST_KIND}" MATCHES "Jali" )

    # In this case, we need the Jali target definition
    if (NOT TARGET Jali )
      message(FATAL_ERROR "Can not define an Jali test before defining Jali binary")
    endif()  

    get_target_property(base Jali OUTPUT_NAME)
    get_target_property(dir  Jali OUTPUT_DIRECTORY)
    #set(test_exec "${dir}/${base}")
    set(test_exec "${SSC_BINARY_DIR}/${base}")
   
  else() 
   
    # Do not set if this variable is empty
    if ( Jali_TEST_UNPARSED_ARGUMENTS )
      list(GET Jali_TEST_UNPARSED_ARGUMENTS 0 test_exec)
      list(REMOVE_AT Jali_TEST_UNPARSED_ARGUMENTS 0)
    endif()  

    # Throw an error if test_exec is not defined
    if ( NOT test_exec )
      message(FATAL_ERROR "Must define a test executable for ${test_name}")
    endif()  

    # Create the executable if SOURCE is defined
    if(Jali_TEST_SOURCE)
      add_executable(${test_exec} ${Jali_TEST_SOURCE})
      set_target_properties(${test_exec} PROPERTIES COMPILE_FLAGS "-DCMAKE_BINARY_DIR=\\\"${CMAKE_BINARY_DIR}\\\" -DCMAKE_SOURCE_DIR=\\\"${CMAKE_SOURCE_DIR}\\\"")
    endif()

    # Add link libraries if needed
    if(Jali_TEST_LINK_LIBS)
      target_link_libraries(${test_exec} ${Jali_TEST_LINK_LIBS})
    endif()

  endif()  

  
  # --- Define the test arguments

  set(test_args "${Jali_TEST_UNPARSED_ARGUMENTS}")
  if ( "${Jali_TEST_KIND}" MATCHES "Jali" )
    
    # In this case, we need an Jali input file
    if ( NOT Jali_TEST_Jali_INPUT )
      message(FATAL_ERROR "Jali tests require an Jali input file")
    endif()

    set(test_args "--xml_file=${Jali_TEST_Jali_INPUT};${test_args}")

  endif()  

  # --- Add test

  # Adjust the execuable name if NOT fullpath AND TESTS_REQUIRE_FULLPATH is set
  if ( TESTS_REQUIRE_FULLPATH )
    if ( NOT ("${test_exec}" MATCHES "^/") )
      set(_tmp      "${CMAKE_CURRENT_BINARY_DIR}/${test_exec}")
      set(test_exec "${_tmp}")
    endif()  
  endif()

  # Construct the test execution command
  set(add_test_exec)
  set(add_test_args)
  if (Jali_TEST_PARALLEL)

    if ( MPI_EXEC_GLOBAL_ARGS )
      separate_arguments(global_mpi_args UNIX_COMMAND "${MPI_EXEC_GLOBAL_ARGS}")
    endif() 

    set(add_test_exec ${MPI_EXEC})
    set(add_test_args
                      ${MPI_EXEC_NUMPROCS_FLAG}
                      ${Jali_TEST_NPROCS}
                      ${global_mpi_args}
                      ${Jali_TEST_MPI_EXEC_ARGS}
                      ${MPI_EXEC_PREFLAGS}
                      ${test_exec}
                      ${MPI_EXEC_POSTFLAGS}
                      ${test_args})
  else()
    set(add_test_exec ${test_exec})
    set(add_test_args ${test_args})
  endif()

  # Call add_test
  add_test(NAME ${test_name} COMMAND ${add_test_exec} ${add_test_args})

  # --- Add test properties

  # Labels
  _add_test_kind_label(${test_name} ${Jali_TEST_KIND})
  if ( Jali_TEST_PARALLEL AND Jali_TEST_NPROCS )
    if ( ${Jali_TEST_NPROCS} GREATER 1 )
      _append_test_label(${test_name} PARALLEL)
    else()
      _append_test_label(${test_name} SERIAL)
    endif()  
  else()  
    _append_test_label(${test_name} SERIAL)
  endif()

  # Add test dependencies
  if ( Jali_TEST_DEPENDS )
    set_tests_properties(${test_name} PROPERTIES
                         DEPENDS "${Jali_TEST_DEPENDS}")
  endif()		       
  
  # Remaining properties are single valued. Building 
  # test_properties as a list should get past the CMake parser.

  # Timeout
  if ( TESTS_TIMEOUT_THRESHOLD )
    list(APPEND test_properties TIMEOUT ${TESTS_TIMEOUT_THRESHOLD})
  endif()

  # CTest needs to know how procs this test needs
  if ( Jali_TEST_PARALLEL )
    list(APPEND test_properties PROCESSORS ${Jali_TEST_NPROCS})
  endif()

  # Set expected failure flag
  if ( Jali_TEST_EXPECTED_FAIL )
     list(APPEND test_properties WILL_FAIL TRUE)
  endif() 

  if ( test_properties )
    set_tests_properties(${test_name} PROPERTIES ${test_properties})
  endif()

endfunction(ADD_Jali_TEST)




