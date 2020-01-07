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
#    Mesh 
#
#

# Jali module, include files found in JALI_MODULE_PATH
# include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   GEOMETRY_SOURCE_DIR
#   GEOMETRY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOMETRY)


set(GEOMETRY_headers
  BoxRegion.hh
  Domain.hh
  GeometricModel.hh
  Geometry.hh
  LabeledSetRegion.hh
  LogicalRegion.hh
  PlaneRegion.hh
  Point.hh
  PointRegion.hh
  PolygonRegion.hh
  Region.hh
  )
list(TRANSFORM GEOMETRY_headers PREPEND "${GEOMETRY_SOURCE_DIR}/")

set(GEOMETRY_sources
  BoxRegion.cc
  Domain.cc
  GeometricModel.cc
  Geometry.cc
  LabeledSetRegion.cc
  LogicalRegion.cc
  PlaneRegion.cc
  PointRegion.cc
  Region.cc
  )

#
# Library Target: geometry
#
add_library(jali_geometry ${GEOMETRY_sources})
set_target_properties(jali_geometry PROPERTIES PUBLIC_HEADER "${GEOMETRY_headers}")

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(Jali::jali_geometry ALIAS jali_geometry)

# We should change install_interface to include/error_handling but that requires
# changing all client codes
target_include_directories(jali_geometry PUBLIC
  $<BUILD_INTERFACE:${GEOMETRY_BINARY_DIR}>
  $<BUILD_INTERFACE:${GEOMETRY_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Make the error handling library a dependency of this target
target_link_libraries(jali_geometry PUBLIC jali_error_handling)

#
# Unit tests
#
if (BUILD_TESTS)
    
  target_include_directories(jali_geometry PRIVATE ${UnitTest_INCLUDE_DIRS})
  target_link_libraries(jali_geometry PRIVATE ${UnitTest_LIBRARIES})

  # Need to copy files for the tests 
  if (NOT (${GEOMETRY_SOURCE_DIR} EQUAL ${GEOMETRY_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${GEOMETRY_SOURCE_DIR}/test ${GEOMETRY_BINARY_DIR}/test) 
   endif()

   # Test: test point creation and deletion
   add_Jali_test(point-tests test_points
     KIND unit
     SOURCE test/Main.cc test/test_points.cc
     LINK_LIBS jali_geometry ${UnitTest++_LIBRARIES})
   
   # Test: test geometric operators
   add_Jali_test(geometry-operators test_geomops
     KIND unit
     SOURCE test/Main.cc test/test_geometric_ops.cc
     LINK_LIBS jali_geometry ${UnitTest++_LIBRARIES})

   # Test: test region creation - Why aren't we testing these?? Its
   # possible if we switch from using the XML parser to something else
   # to read the specification or just hard code the parameters.
   
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
#                          jali_geometry
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
#                          jali_geometry
#                          ${UnitTest_LIBRARIES}
#                          ${Teuchos_LIBRARIES}
#                          ${Epetra_LIBRARIES})  

endif()

