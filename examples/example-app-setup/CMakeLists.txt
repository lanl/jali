project (JALIAPP)

cmake_minimum_required(VERSION 2.8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#
# Jali Libraries and includes
#

if (NOT Jali_DIR)
  message(FATAL_ERROR "Error: Jali top level installation dir must be defined")
endif()

# Look for the Jali package

find_package(Jali)
if (Jali_FOUND) 
   message(STATUS "Located Jali")
   message(STATUS "Jali_DIR=${Jali_DIR}")
   message(STATUS "Jali_INCLUDE_DIRS=${Jali_INCLUDE_DIRS}")
   message(STATUS "Jali_LIBRARY_DIRS=${Jali_LIBRARY_DIRS}")
   message(STATUS "Jali_LIBRARIES=${Jali_LIBRARIES}")
   message(STATUS "Jali_TPL_INCLUDE_DIRS=${Jali_TPL_INCLUDE_DIRS}")
   message(STATUS "Jali_TPL_LIBRARY_DIRS=${Jali_TPL_LIBRARY_DIRS}")
   message(STATUS "Jali_TPL_LIBRARIES=${Jali_TPL_LIBRARIES}")
endif()


# Default build type is Release (Optimized)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()


# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${PROJECT_SOURCE_DIR})
endif ()

# if its an older CMAKE we have to explicitly enable C++11 support
if ("${CMAKE_MAJOR_VERSION}" LESS 3 OR "{CMAKE_PATCH_VERSION}" LESS 0)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif()


# include files and source files  
list(REVERSE Jali_LIBRARIES)

foreach (Jali_library ${Jali_LIBRARIES})
  message(STATUS "Adding Jali_library ${Jali_library} as an imported library")
  add_library(${Jali_library} STATIC IMPORTED)
  set_property(TARGET ${Jali_library} PROPERTY IMPORTED_LOCATION ${Jali_LIBRARY_DIRS}/lib${Jali_library}.a)
endforeach()

add_executable(testmesh src/main.cc)
include_directories(${PROJECT_SOURCE_DIR}/include ${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})
target_link_libraries(testmesh ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

if ("${CMAKE_MAJOR_VERSION}" GREATER 2 AND "{CMAKE_PATCH_VERSION}" GREATER 0)
  set_property(TARGET testmesh PROPERTY CXX_STANDARD 11)
  set_property(TARGET testmesh PROPERTY CXX_STANDARD_REQUIRED ON)
endif()

