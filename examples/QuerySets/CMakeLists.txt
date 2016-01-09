# -*- mode: cmake -*-

#
#  jali
#    examples
#      QuerySets
#

project(QuerySets)

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

if (NOT (${QuerySets_SOURCE_DIR} EQUAL ${QuerySets_BINARY_DIR}))
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy 
    ${QuerySets_SOURCE_DIR}/hex_3x3x3_sets.exo ${QuerySets_BINARY_DIR}/hex_3x3x3_sets.exo)
endif()

add_executable(QuerySets QuerySets.cc)
target_link_libraries(QuerySets geometry mesh mstk_mesh mesh_factory)

  
