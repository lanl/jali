# -*- mode: cmake -*-

#
#  jali
#    examples
#      DIffusion1D
#

# Jali include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(Diffusion1D Diffusion1D.cc)
target_link_libraries(Diffusion1D mesh mstk_mesh mesh_factory)

  
