# -*- mode: cmake -*-
# 
# Jali Third Party Library (TPL) Definitions
#

# Use <package_name>_ROOT variables in find_package (introduced in 3.12)
if (${CMAKE_VERSION} VERSION_GREATER 3.12)
  cmake_policy(SET CMP0074 NEW)
endif(${CMAKE_VERSION} VERSION_GREATER 3.12)

# Standard CMake modules see CMAKE_ROOT/Modules
include(FeatureSummary)

# Jali CMake modules see <root source>/tools/cmake
include(AddImportedLibrary)
include(TrilinosMacros)
include(PrintVariable)

##############################################################################
# ------------------------ Required Libraries -------------------------------#
##############################################################################

##############################################################################
# MPI
##############################################################################
include(CheckMPISourceCompiles)
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)

if ( NOT MPI_WRAPPERS_IN_USE )
  find_package(MPI REQUIRED)
  
  # Warn the user if MPI information is not found
  if (NOT MPI_C_FOUND)
    message(WARNING "Failed to locate MPI C include and library files")
  endif()
  if (NOT MPI_CXX_FOUND)
    message(WARNING "Failed to locate MPI C++ include and library files")
  endif()

  # We have to compile Jali and Jali TPLs with MPI - so forcibly set the compilers to the MPI wrappers 
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "MPI C wrapper to use" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "MPI C++ wrapper to use" FORCE)
  
endif ( NOT MPI_WRAPPERS_IN_USE )

# - Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if (NOT MPI_EXEC_NUMPROCS_FLAG)
  if (MPIEXEC_NUMPROC_FLAG)
    set(MPI_EXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "Set MPI number of procs flag from FindMPI")
  else()
    set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
  endif()
endif()
  
# - Maximum number of processors. This is a limit for the test suite
#   Some tests require too many processors and it increases the execution time
#   considerably. 
set(MPI_EXEC_MAX_NUMPROCS_DFLT 8)
if (NOT MPI_EXEC_MAX_NUMPROCS)
  include(ProcessorCount)
  ProcessorCount(proc_count)
  if (NOT proc_count EQUAL 0)
    math(EXPR MPI_EXEC_MAX_NUMPROCS "${proc_count} * 2") 
    message(STATUS "Detected ${proc_count} processors and will set maximum to ${MPI_EXEC_MAX_NUMPROCS}")
  else()
    set(MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DFLT})
  endif()
endif()  

# - Set the pre and post flags
#   Usage:
#   ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} PROCS ${MPI_EXEC_PREFLAGS} EXECUTABLE ${MPI_EXEC_POSTFLAGS}
if (NOT MPI_EXEC_PREFLAGS)
  if (MPIEXEC_PREFLAGS)
    set(MPI_EXEC_PREFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute pre flags")
  endif()
endif()

if (NOT MPI_EXEC_POSTFLAGS)
  if (MPIEXEC_POSTFLAGS)
    set(MPI_EXEC_POSTFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute post flags")
  endif()
endif()

##############################################################################
# Boost
##############################################################################

# CMake 2.8.6 FindBoost stops at version 1.46
# Add more versions to the search see cmake --help-module FindBoost for
# more information.
set(Boost_ADDITIONAL_VERSIONS 
    1.47 1.47.0
    1.48 1.48.0
    1.49 1.49.0
    1.50 1.50.0
    1.51 1.51.0
    1.52 1.52.0
    1.53 1.53.0
    1.54 1.55.0)
find_package( Boost COMPONENTS system filesystem program_options regex REQUIRED)
option(ENABLE_BOOST "Boost info" ON)
add_feature_info(Boost ENABLE_BOOST "Cpp Extension library http://www.boost.org")

if ( Boost_VERSION) 

  if ( ${Boost_VERSION} VERSION_LESS 1.63 )
    message(WARNING "Found Boost version ${Boost_VERSION} which"
                    " is older than the supported (1.46) version.")
  endif()

  # The Boost filesystem library changed and deprecated some functions.
  # This define should be used when packages include boost/filesystem.hpp
  # and packages any of these new or deprecated functions.
  # The change from version 2 to 3 occurred with the 1.49 Boost release.
  # Please refer to the online documentation at www.boost.org.
  if ( "${Boost_VERSION}" VERSION_LESS "1.34" )
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=1")
  elseif ( "${Boost_VERSION}" VERSION_LESS "1.49" )
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=2")
  else()
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=3")
  endif()  
endif()

##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################

find_package(HDF5 1.8.18 REQUIRED C HL)
if ( NOT HDF5_IS_PARALLEL ) 
    message(WARNING     "The HDF5 installation found in ${HDF5_DIR} is not "
                        "a parallel build. At this time, this installation "
                        "is compatible with other TPLs. Soon Jali will "
                        "require a parallel enabled HDF5. Please update your "
                        "HDF5 installation to include MPI I/O symbols"
            )            
endif(NOT HDF5_IS_PARALLEL)
option(ENABLE_HDF5 "HDF5 info" ON)
add_feature_info(HDF5 ENABLE_HDF5 "I/O library that creates HDF5 formatted files http://www.hdfgroup.org/HDF5")


##############################################################################
# Trilinos http://trilinos.sandia.gov
##############################################################################
# This command alters Trilinos_DIR. If it finds the configuration file
# Trilinos_DIR is set to the path the configuration file was found.
if ( NOT Trilinos_INSTALL_PREFIX )
  message(WARNING "Use Trilinos_INSTALL_PREFIX"
                  " to define the Trilinos installation location"
		  "\n-DTrilinos_INSTALL_PREFIX:PATH=<trilnos directory>\n")
endif()
set(Trilinos_MINIMUM_VERSION 12.0.0)
find_package(Trilinos ${Trilinos_MINIMUM_VERSION} REQUIRED
             PATHS ${Trilinos_INSTALL_PREFIX}
             PATH_SUFFIXES lib/cmake/Trilinos)
            
if (Trilinos_FOUND)
  message(STATUS "Found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
  trilinos_package_enabled_tpls(Trilinos)           

  if ("${Trilinos_VERSION}" VERSION_LESS ${Trilinos_MINIMUM_VERSION}) 
    message(FATAL_ERROR "Trilinos version ${Trilinos_VERSION} is not sufficient."
                        " Amanzi requires at least version ${Trilinos_MINIMUM_VERSION}")
  endif()

    message(STATUS "Found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
    trilinos_package_enabled_tpls(Trilinos)           

    if ( "${Trilinos_VERSION}" VERSION_LESS ${Trilinos_MINIMUM_VERSION} ) 
      message(FATAL_ERROR "Trilinos version ${Trilinos_VERSION} is not sufficient."
	                  " Jali requires at least version ${Trilinos_MINIMUM_VERSION}")
    endif()

    # Now check optional Trilinos packages

    # STK - mesh 
    if ( ENABLE_STK_Mesh )
      find_package(STK
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib)
      if (STK_FOUND)
	message(STATUS "Located Trilinos package STK: ${STK_DIR}")
        trilinos_package_enabled_tpls(STK)
	foreach( _inc "${STK_TPL_INCLUDE_DIRS}")
	    list(APPEND STK_INCLUDE_DIRS "${_inc}")
	endforeach()
      else()  
        message(WARNING "Could not locate STK in ${Trilinos_DIR}. Will not enable STK_Mesh")
        set(ENABLE_STK_Mesh OFF CACHE BOOL "Disable STK Mesh capability" FORCE)
      endif()

    endif()

    # Zoltan - required by MSTK mesh class 
    if ( ENABLE_MSTK_Mesh )
      find_package(Zoltan
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib)
      if (Zoltan_FOUND)
	message(STATUS "Located Trilinos package Zoltan: ${Zoltan_DIR}")
        trilinos_package_enabled_tpls(Zoltan)
	foreach( _inc "${ZOLTAN_TPL_INCLUDE_DIRS}")
	    list(APPEND ZOLTAN_INCLUDE_DIRS "${_inc}")
	endforeach()
      else()  
        message(WARNING "Could not locate Zoltan in ${Trilinos_DIR}. Will not enable MSTK_Mesh")
        set(ENABLE_MSTK_Mesh OFF CACHE BOOL "Disable MSTK Mesh capability" FORCE)
      endif()

    endif()


    # Now update the Trilinos_LIBRARIES and INCLUDE_DIRS
    foreach( _inc "${Trilinos_TPL_INCLUDE_DIRS}")
      list(APPEND Trilinos_INCLUDE_DIRS "${_inc}")
      list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
    endforeach()

else()
    message(FATAL_ERROR "Can not locate Trilinos configuration file\n"
                        " Please define the location of your Trilinos installation\n"
                        "using -D Trilinos_DIR:FILEPATH=<install path>\n")
endif()    

##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################
find_package(NetCDF REQUIRED)
option(ENABLE_NetCDF "NetCDF info" ON)
add_feature_info(NetCDF ENABLE_NetCDF "Network Common Data Format (NetCDF) http://www.unidata.ucar.edu/software/netcdf/")


##############################################################################
# Exodus II -http://sourceforge.net/projects/exodusii
##############################################################################
find_package(ExodusII REQUIRED)
option(ENABLE_ExodusII "ExodusII info" ON)
add_feature_info(ExodusII ENABLE_ExodusII "File format library from Sandia National Labs. https://github.com/gsjaardema/seacas")




##############################################################################
#---------------------------- Mesh Frameworks -------------------------------#
##############################################################################


##############################################################################
# STK - Sierra Mesh Tool Kit part of Trilinos
##############################################################################
option(ENABLE_STK_Mesh  "Build Jali with the STK mesh framework" OFF)
add_feature_info(STK_Mesh
                 ENABLE_STK_Mesh
                 "Sierra Mesh Tool Kit (STK Mesh) a Trilinos package"
                 )


##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################
option(ENABLE_MOAB_Mesh "Build Jali with the MOAB mesh framework" OFF)
add_feature_info(MOAB_Mesh
                 ENABLE_MOAB_Mesh
                 "A Mesh-Oriented datABase"
                 )
if (ENABLE_MOAB_Mesh)
    find_package(MOAB REQUIRED)
endif()

##############################################################################
# MSTK - https://github.com/MeshToolkit/mstk
##############################################################################
option(ENABLE_MSTK_Mesh "Build Jali with the MOAB mesh framework" OFF)
add_feature_info(MSTK_Mesh
                 ENABLE_MSTK_Mesh
                 "A mesh framework"
                 )
if (ENABLE_MSTK_Mesh)
  find_package(MSTK REQUIRED)
endif() 





##############################################################################
#-------------------------- Optional Libraries ------------------------------#
##############################################################################

##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################
option(ENABLE_UnitTest "Build Jali unit tests. Requires UnitTest++" ON)
add_feature_info(UnitTest
                 ENABLE_UnitTest
                 "C++ unit test framework"
                 )
if (ENABLE_UnitTest)
    find_package(UnitTest)
endif()    

##############################################################################
# OpenMP - http://openmp.org/
#
# comment out add_feature_info per
# https://software.lanl.gov/ascem/trac/ticket/413#comment:1
##############################################################################
option(ENABLE_OpenMP "Build Jali executables with OpenMP" OFF)
#add_feature_info(OpenMP
#                 ENABLE_OpenMP
#                 "OpenMP, multi-platform shared-memory parallel programming"
#                 )
if (ENABLE_OpenMP)
    find_package(OpenMP)
endif()

