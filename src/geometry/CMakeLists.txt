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
                   LINK_LIBS geometry ${UnitTest_LIBRARIES})

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

