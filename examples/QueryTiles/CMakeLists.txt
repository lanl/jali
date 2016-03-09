# -*- mode: cmake -*-

#
#  jali
#    examples
#      QueryTopology
#

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(QueryTiles QueryTiles.cc)
target_link_libraries(QueryTiles mesh mstk_mesh mesh_factory)

  
