# -*- mode: cmake -*-

#
#  Jali
#    Mesh 
#
#

# Jali module, include files found in JALI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   GEOMETRY_SOURCE_DIR
#   GEOMETRY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOMETRY)

# Jali include directories
include_directories(${DBC_SOURCE_DIR} ${FUNCS_SOURCE_DIR} ${ATK_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# geometry library
#
file (GLOB geometry_source_files "*.cc")
file(GLOB  geometry_inc_files "*.h*")
#add_Jali_library(geometry 
#                   SOURCE ${geometry_source_files} HEADERS ${geometry_inc_files}
#                   LINK_LIBS error_handling functions ${Teuchos_LIBRARIES} )

add_Jali_library(geometry 
                   SOURCE ${geometry_source_files} HEADERS ${geometry_inc_files}
                   LINK_LIBS error_handling)


#
# Unit tests
#
if (BUILD_TESTS)
    
   # Add UnitTest ExodusII include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Need to copy files for the tests 
   if (NOT (${GEOMETRY_SOURCE_DIR} EQUAL ${GEOMETRY_BINARY_DIR}))
       execute_process(COMMAND ${CMAKE_COMMAND} -E 
         copy_directory ${GEOMETRY_SOURCE_DIR}/test ${GEOMETRY_BINARY_DIR}/test) 
   endif()

   # Test: test point creation and deletion
   add_Jali_test(point-tests test_points
                   KIND unit
		   SOURCE test/Main.cc test/test_points.cc
		   LINK_LIBS geometry ${UnitTest_LIBRARIES})

   # Test: test geometric operators
   add_Jali_test(geometry-operators test_geomops
                   KIND unit
                   SOURCE test/Main.cc test/test_geometric_ops.cc
                   LINK_LIBS geometry ${Teuchos_LIBRARIES} ${UnitTest_LIBRARIES})

   # Test: test region creation
#   add_Jali_test(geometry-region test_region
#                   KIND unit
#                   SOURCE  
#                          test/Main.cc 
#                          test/test_box_region.cc
#                          test/test_plane_region.cc
#                          test/test_labeledset.cc
#                          test/test_colorfunc_region.cc
#			  test/test_polygon_region.cc
#                   LINK_LIBS 
#                          geometry
#                          ${UnitTest_LIBRARIES}
#                          ${Teuchos_LIBRARIES}
#                          ${Epetra_LIBRARIES})  

   # Test: test geometric model creation
#   add_Jali_test(geometry-model test_geometric_model
#                   KIND unit
#                   SOURCE
#                          test/Main.cc
#                          test/test_geometric_model.cc
#                   LINK_LIBS 
#                          geometry
#                          ${UnitTest_LIBRARIES}
#                          ${Teuchos_LIBRARIES}
#                          ${Epetra_LIBRARIES})  
           

endif()
