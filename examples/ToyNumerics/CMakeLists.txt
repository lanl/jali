# -*- mode: cmake -*-

#
#  jali
#    examples
#      ToyNumerics
#

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${JALI_STATE_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(ToyNumerics ToyNumerics.cc)
target_link_libraries(ToyNumerics mesh mstk_mesh mesh_factory jali_state)

  
