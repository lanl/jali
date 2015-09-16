# -*- mode: cmake -*-

#
#  Jali
#    Mesh 
#
#

# Jali module, include files found in JALI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH)

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh
#
add_Jali_library(mesh SOURCE Mesh.cc LINK_LIBS geometry)

#
# Library: mesh_audit
#
#add_Jali_library(mesh_audit SOURCE MeshAudit.cc LINK_LIBS error_handling)

#
# Install Header files
#
file(GLOB mesh_inc_files "*.hh")
add_install_include_file(${mesh_inc_files})



# Derived and support classes

add_subdirectory(mesh_simple)

# Mesh Frameworks

# STK (Trilinos Package)
#if (ENABLE_STK_Mesh)
#    add_subdirectory(mesh_stk)
#endif(ENABLE_STK_Mesh)    

# MOAB Mesh framework
#if (ENABLE_MOAB_Mesh)
#    add_subdirectory(mesh_moab)
#endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
    add_subdirectory(mesh_mstk)
endif(ENABLE_MSTK_Mesh)    


# Factory class

add_subdirectory(mesh_factory)


if (BUILD_TESTS)

    # Include dirs for test
    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})

    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_SOURCE_DIR} EQUAL ${MESH_BINARY_DIR}))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${MESH_SOURCE_DIR}/test ${MESH_BINARY_DIR}/test) 
    endif()

    # Test link libraries
    set(test_link_libs geometry simple_mesh mesh_factory ${UnitTest_LIBRARIES})

    # STK (Trilinos Package)
   if (ENABLE_STK_Mesh)
     list(APPEND test_link_libs stk_mesh)
   endif(ENABLE_STK_Mesh)    

    # MOAB Mesh framework
   if (ENABLE_MOAB_Mesh)
     list(APPEND test_link_libs moab_mesh)
   endif(ENABLE_MOAB_Mesh)    
    
    # MSTK Mesh framework
    if (ENABLE_MSTK_Mesh)
      list(APPEND test_link_libs mstk_mesh)
    endif(ENABLE_MSTK_Mesh)    

    # Test: mesh geometry

    set(test_src_files test/Main.cc test/test_mesh_geometry.cc)

    add_Jali_test(mesh_geometry test_mesh_geometry
                    KIND unit
		    SOURCE ${test_src_files}
		    LINK_LIBS ${test_link_libs})

    # Test mesh entity iterators

    add_Jali_test(mesh_iterators test_entity_iterators
                  KIND unit
                  SOURCE test/Main.cc test/test_entity_iterators.cc
		  LINK_LIBS ${test_link_libs})
endif()
  
