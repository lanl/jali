# -*- mode: cmake -*-

#
#  Jali
#    MSTK Mesh frameworks
#

# Jali module, include files found in JALI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MSTK_SOURCE_DIR
#   MSTK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MSTK)

# External (TPL) include directories

# Reverse the items in MSTK_INCLUDE_DIRS because when we ask the
# include_directories command to prepend a set of directories using
# the BEFORE keyword, it prepends each item in the MSTK_INCLUDE_DIRS
# list, thereby effectively reversing the order of the
# directories. For example, if we had MYDIRS=a;b;c and we specified
# MYDIRS in the include_directories command, then the
# INCLUDE_DIRECTORIES property for will end up with
# c;b;a;whatever_existed_before

list(REVERSE MSTK_INCLUDE_DIRS)
include_directories(BEFORE ${MSTK_INCLUDE_DIRS})

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})


# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

#
# Library: mstk_mesh
#
set(mstk_mesh_files Mesh_MSTK.cc) 
file(GLOB mstk_inc_files "*.hh")
add_Jali_library(mstk_mesh
                   SOURCE ${mstk_mesh_files} HEADERS ${mstk_inc_files}
                   LINK_LIBS geometry mesh ${MSTK_LIBRARIES} ${Zoltan_LIBRARIES})

if (BUILD_TESTS)
    
    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MSTK_SOURCE_DIR}/test ${MSTK_BINARY_DIR}/test) 
    endif()


    # Test: mstk_mesh
    add_Jali_test(mstk_mesh test_mstk_mesh
                    KIND int
                    SOURCE
                         test/Main.cc 
                         test/test_hex_1x1x1.cc
                         test/test_hex_2x2x1.cc
                         test/test_hex_3x3x3.cc
                         test/test_hex_gen_3x3x3.cc
                         test/test_quad_gen_3x3.cc
			 test/test_write_read_fields.cc
			 test/test_block_partition.cc
                    LINK_LIBS mstk_mesh ${UnitTest_LIBRARIES}) 

    # Test: mstk_mesh_parallel
    add_Jali_test(mstk_mesh_parallel test_parallel_mstk_mesh
                    KIND unit
                    NPROCS 4
                    SOURCE
                         test/Main.cc
                         test/test_hex_3x3x3_4P.cc
    			 test/test_hex_3x3x3_par_read_4P.cc
                         test/test_quad_gen_3x3_4P.cc
                         test/test_hex_gen_3x3x3_4P.cc
			 test/test_edges_4P.cc
                    LINK_LIBS mstk_mesh ${UnitTest_LIBRARIES})

endif()

