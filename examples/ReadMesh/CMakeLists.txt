# -*- mode: cmake -*-

#
#  jali
#    examples
#      ReadMesh
#

project(ReadMesh)

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

if (NOT (${ReadMesh_SOURCE_DIR} EQUAL ${ReadMesh_BINARY_DIR}))
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy 
    ${ReadMesh_SOURCE_DIR}/quadtri.exo ${ReadMesh_BINARY_DIR}/quadtri.exo)
endif()

add_executable(ReadMesh ReadMesh.cc)
target_link_libraries(ReadMesh mesh mstk_mesh mesh_factory)

  
