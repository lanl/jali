# -*- mode: cmake -*-

#
#  jali
#    statemanager
#

# Jali module, include files found in JALI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   JALI_STATE_SOURCE_DIR
#   JALI_STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(JALI_STATE)

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

#
# Library: 
#
add_Jali_library(jali_state 
  SOURCE JaliStateVector.h JaliStateVector.cc JaliState.h JaliState.cc 
  LINK_LIBS mesh)

#
# Install Header files
#
file(GLOB state_inc_files "*.h")
add_install_include_file(${state_inc_files})


if (BUILD_TESTS)

    # Include dirs for test
    include_directories(${JALI_STATE_SOURCE_DIR})
    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})

    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${JALI_STATE_SOURCE_DIR} EQUAL ${JALI_STATE_BINARY_DIR}))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${JALI_STATE_SOURCE_DIR}/test ${JALI_STATE_BINARY_DIR}/test) 
    endif()

    # Test link libraries
    set(test_link_libs jali_state mesh mstk_mesh mesh_factory ${UnitTest_LIBRARIES})

    # Test state vectors

    set(test_src_files test/Main.cc test/test_jali_state_vector.cc)

    add_Jali_test(jali_state_vectors test_jali_state_vectors
                    KIND unit
		    SOURCE ${test_src_files}
		    LINK_LIBS ${test_link_libs})

    # Test state

    set(test_src_files test/Main.cc test/test_jali_state.cc)

    add_Jali_test(jali_state test_jali_state
                  KIND unit
		  SOURCE ${test_src_files}
		  LINK_LIBS ${test_link_libs})

endif()
  
