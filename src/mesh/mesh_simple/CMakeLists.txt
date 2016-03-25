# -*- mode: cmake -*-

#
#  Jali
#    Mesh Base Class
#

# Jali module, include files found in JALI_MODULE_PATH
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SIMPLE_SOURCE_DIR
#   MESH_SIMPLE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_SIMPLE)

# Jali include directories
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories

# Library: simple_mesh
file(GLOB simple_mesh_source_files "*.cc")
file(GLOB simple_inc_files "*.hh")
add_Jali_library(simple_mesh
                   SOURCE ${simple_mesh_source_files} HEADERS ${simple_inc_files}
                   LINK_LIBS mesh geometry error_handling)

if (BUILD_TESTS)

    # Add UnitTest include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test uses Epetra

   # Need to copy files for the tests 
    if (NOT (${MESH_SIMPLE_SOURCE_DIR} EQUAL ${MESH_SIMPLE_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MESH_SIMPLE_SOURCE_DIR}/test ${MESH_SIMPLE_BINARY_DIR}/test) 
    endif()

    # Test: simple_mesh
    add_Jali_test(simple_mesh test_simple_mesh
                    KIND unit
		    SOURCE
                          test/Main.cc 
                          test/test_maps.cc 
                          test/test_cell_numbering.cc 
                          test/test_face_adj_cells.cc 
                          test/test_node_cell_faces.cc
			  test/test_geometry.cc
                    LINK_LIBS simple_mesh ${UnitTest_LIBRARIES})

endif()

