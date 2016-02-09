# -*- mode: cmake -*-
# 
# Jali Third Party Library (TPL) Definitions
#

# Standard CMake modules see CMAKE_ROOT/Modules
include(FeatureSummary)

# Jali CMake modules see <root source>/tools/cmake
include(CheckMPISourceCompiles)
include(TrilinosMacros)
include(PrintVariable)


##############################################################################
# ------------------------ Required Libraries -------------------------------#
##############################################################################

##############################################################################
# MPI
##############################################################################
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)

if ( NOT MPI_WRAPPERS_IN_USE )
  message(WARNING "At this time, Jali must be compiled with MPI wrappers."
                  " Build will likely fail. Please define CMAKE_*_COMPILER"
		  " parameters as MPI compiler wrappers and re-run cmake.")
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
set_feature_info(Boost
                 "C++ Extension library"
                 "http://www.boost.org"
                 "Required by the MPC")

if ( Boost_VERSION) 

  if ( ${Boost_VERSION} VERSION_LESS 1.46 )
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

# We need to use the project-local HDF5 finder. Temporarily Change
# policy CMP0017 if were using cmake 2.8.3 or later
if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 OLD)
endif()

find_package(HDF5 1.8.0 REQUIRED)
if ( NOT HDF5_IS_PARALLEL ) 
    message(WARNING     "The HDF5 installation found in ${HDF5_DIR} is not "
                        "a parallel build. At this time, this installation "
                        "is compatible with other TPLs. Soon Jali will "
                        "require a parallel enabled HDF5. Please update your "
                        "HDF5 installation to include MPI I/O symbols"
            )            
endif(NOT HDF5_IS_PARALLEL)
set_feature_info(HDF5
                "I/O library that creates HDF5 formatted files"
                "http://www.hdfgroup.org/HDF5"
                "Required library for several components in Jali"
                )

# Restore policy of preferring offical CMake modules over local ones.
if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

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
set(Trilinos_MINIMUM_VERSION 11.0.3)
find_package(Trilinos 
             PATHS ${Trilinos_INSTALL_PREFIX}
             PATH_SUFFIXES include)
            
if ( Trilinos_FOUND )

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
set_feature_info(NetCDF
                 "Network Common Data Format (NetCDF)"
                 "http://www.unidata.ucar.edu/software/netcdf/"
                 "Required by ExodusII library")


##############################################################################
# Exodus II -http://sourceforge.net/projects/exodusii
##############################################################################
find_package(ExodusII REQUIRED)
set_feature_info(ExodusII
                 "File format library. Originated from Sandia."
                 "http://sourceforge.net/projects/exodusii/"
                 "Required by all the mesh frameworks to read mesh files")


##############################################################################
# XERCES-C - http://http://xerces.apache.org/xerces-c/
##############################################################################
#find_package(XERCES REQUIRED)
#set_feature_info(XERCES
#	         "Validating XML parser")


##############################################################################
############################ Option Processing ###############################
##############################################################################







##############################################################################
#---------------------------- Mesh Frameworks -------------------------------#
##############################################################################

# Enable ALL possible mesh frameworks
#option(ENABLE_ALL_Mesh "Build all Jali mesh frameworks" OFF)
#if(ENABLE_ALL_Mesh)
#    set(ENABLE_STK_Mesh ON)
#    set(ENABLE_MOAB_Mesh ON)
#    set(ENABLE_MSTK_Mesh ON)
#endif()    
#set_feature_info(ALL_Mesh
#                 ENABLE_ALL_Mesh
#                 "Build all available mesh frameworks"
#                  )    

##############################################################################
# STK - Sierra Mesh Tool Kit part of Trilinos
##############################################################################
option(ENABLE_STK_Mesh  "Build Jali with the STK mesh framework" OFF)
set_feature_info(STK_Mesh
                 ENABLE_STK_Mesh
                 "Sierra Mesh Tool Kit (STK Mesh) a Trilinos package"
                 )


##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################
option(ENABLE_MOAB_Mesh "Build Jali with the MOAB mesh framework" OFF)
set_feature_info(MOAB_Mesh
                 ENABLE_MOAB_Mesh
                 "A Mesh-Oriented datABase"
                 )
if (ENABLE_MOAB_Mesh)
    find_package(MOAB REQUIRED)
endif()

##############################################################################
# MSTK - https://software.lanl.gov/MeshTools/trac/raw-attachment/wiki/WikiStart/mstk-1.80.tar.gz
##############################################################################
option(ENABLE_MSTK_Mesh "Build Jali with the MOAB mesh framework" OFF)
set_feature_info(MSTK_Mesh
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
set_feature_info(UnitTest
                 ENABLE_UnitTest
                 "C++ unit test framework"
                 )
if (ENABLE_UnitTest)
    find_package(UnitTest)
endif()    

##############################################################################
# OpenMP - http://openmp.org/
#
# comment out set_feature_info per
# https://software.lanl.gov/ascem/trac/ticket/413#comment:1
##############################################################################
option(ENABLE_OpenMP "Build Jali executables with OpenMP" OFF)
#set_feature_info(OpenMP
#                 ENABLE_OpenMP
#                 "OpenMP, multi-platform shared-memory parallel programming"
#                 )
if (ENABLE_OpenMP)
    find_package(OpenMP)
    find_package(OpenMP_Fortran)
endif()

