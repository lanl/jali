#  -*- mode: cmake -*-

#
# TPLVersions
#
#    Define the versions, approved download locations for each TPL
#

#
# TPL: Jali Collection of TPLs
#
#   Define a "version number" for the collection of TPLs listed here.
#   It's not clear this is the best way to include this information, 
#   but it's a reasonable place to start.
#   
#   Upgrade History:
#
#   1.0.0       - first version reference used in installations

include(CMakeParseArguments)

MACRO(LIST_LENGTH var)
  SET(entries)
  FOREACH(e ${ARGN})
    SET(entries "${entries}.")
  ENDFOREACH(e)
  STRING(LENGTH "${entries}" ${var})
ENDMACRO(LIST_LENGTH)

# this macro appends version number defines to the tpl_versions.h include file
macro(Jali_tpl_version_write)
  set(singleValueArgs FILENAME PREFIX)
  set(multiValueArgs VERSION)
  set(options "")
  
  cmake_parse_arguments(LOCAL "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  list_length(length ${LOCAL_VERSION})

  if (length GREATER 0) 
    list(GET LOCAL_VERSION 0 MAJOR)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MAJOR ${MAJOR}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MAJOR\n")
  endif()

  if (length GREATER 1)
    list(GET LOCAL_VERSION 1 MINOR)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MINOR ${MINOR}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MINOR\n")
  endif()

  if (length GREATER 2)
    list(GET LOCAL_VERSION 2 PATCH)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_PATCH ${PATCH}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_PATCH\n")
  endif()

  file(APPEND ${LOCAL_FILENAME} "\n")

endmacro(Jali_tpl_version_write)


set (JALI_TPLS_VERSION_MAJOR 1)
set (JALI_TPLS_VERSION_MINOR 0)
set (JALI_TPLS_VERSION_PATCH 0)
set (JALI_TPLS_VERSION ${JALI_TPLS_VERSION}.${JALI_TPLS_VERSION_MINOR}.${JALI_TPLS_VERSION_PATCH})
#   Not sure how to create a meaningful hash key for the collection

#
# TPL: Xerces
#
set(XERCES_VERSION_MAJOR 3)
set(XERCES_VERSION_MINOR 1)
set(XERCES_VERSION_PATCH 1)
set(XERCES_VERSION ${XERCES_VERSION_MAJOR}.${XERCES_VERSION_MINOR}.${XERCES_VERSION_PATCH})
set(XERCES_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(XERCES_ARCHIVE_FILE   xerces-c-${XERCES_VERSION}.tar.gz)
set(XERCES_MD5_SUM        6a8ec45d83c8cfb1584c5a5345cb51ae ) 

#
# TPL: OpenMPI
#
set(OpenMPI_VERSION_MAJOR 1)
set(OpenMPI_VERSION_MINOR 4)
set(OpenMPI_VERSION_PATCH 4)
set(OpenMPI_VERSION ${OpenMPI_VERSION_MAJOR}.${OpenMPI_VERSION_MINOR}.${OpenMPI_VERSION_PATCH})
set(OpenMPI_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(OpenMPI_ARCHIVE_FILE   openmpi-${OpenMPI_VERSION}.tar.bz2)
set(OpenMPI_MD5_SUM        e58a1ea7b8af62453aaa0ddaee5f26a0) 

#
# TPL: CURL
#
set(CURL_VERSION_MAJOR 7)
set(CURL_VERSION_MINOR 37)
set(CURL_VERSION_PATCH 0)
set(CURL_VERSION ${CURL_VERSION_MAJOR}.${CURL_VERSION_MINOR}.${CURL_VERSION_PATCH})
set(CURL_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION}.tar.bz2)
set(CURL_MD5_SUM        7dda0cc2e4136f78d5801ac347be696b)

#
# TPL: zlib
#
set(ZLIB_VERSION_MAJOR 1)
set(ZLIB_VERSION_MINOR 2)
set(ZLIB_VERSION_PATCH 6)
set(ZLIB_VERSION ${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH})
set(ZLIB_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ZLIB_ARCHIVE_FILE   zlib-${ZLIB_VERSION}.tar.gz)
set(ZLIB_MD5_SUM        618e944d7c7cd6521551e30b32322f4a) 

#
# TPL: METIS
#
set(METIS_VERSION_MAJOR 5)
set(METIS_VERSION_MINOR 1)
set(METIS_VERSION_PATCH 0)
set(METIS_VERSION ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
set(METIS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(METIS_ARCHIVE_FILE   metis-${METIS_VERSION}.tar.gz)
set(METIS_MD5_SUM        5465e67079419a69e0116de24fce58fe)

#
# TPL: UnitTest
#
set(UnitTest_VERSION_MAJOR 1)
set(UnitTest_VERSION_MINOR 4)
set(UnitTest_VERSION ${UnitTest_VERSION_MAJOR}.${UnitTest_VERSION_MINOR})
set(UnitTest_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(UnitTest_ARCHIVE_FILE   unittest-cpp-${UnitTest_VERSION}.zip)
set(UnitTest_MD5_SUM       bd373a53403ed51ea1bbb60b1952d7e3) 

#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 58)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.bz2)
set(Boost_MD5_SUM        b8839650e61e9c1c0a89f371dd475546)

#
# TPL: BoostCmake
#
set(BoostCmake_VERSION_MAJOR 1)
set(BoostCmake_VERSION_MINOR 46)
set(BoostCmake_VERSION_PATCH 1)
set(BoostCmake_VERSION        ${BoostCmake_VERSION_MAJOR}.${BoostCmake_VERSION_MINOR}.${BoostCmake_VERSION_PATCH})
set(BoostCmake_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(BoostCmake_ARCHIVE_FILE   boost-cmake-cmake-${BoostCmake_VERSION}.tar.gz)
set(BoostCmake_MD5_SUM        ) 

#
# TPL: HDF5
#
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 8)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(HDF5_ARCHIVE_FILE   hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_MD5_SUM        1196e668f5592bfb50d1de162eb16cff)      

#
# TPL: NetCDF
#
set(NetCDF_VERSION_MAJOR 4)
set(NetCDF_VERSION_MINOR 3)
set(NetCDF_VERSION_PATCH 2)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(NetCDF_ARCHIVE_FILE   netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        2fd2365e1fe9685368cd6ab0ada532a0)

#
# TPL: NetCDF Fortran
#
set(NetCDF_Fortran_VERSION_MAJOR 4)
set(NetCDF_Fortran_VERSION_MINOR 2)
set(NetCDF_Fortran_VERSION ${NetCDF_Fortran_VERSION_MAJOR}.${NetCDF_Fortran_VERSION_MINOR})
set(NetCDF_Fortran_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(NetCDF_Fortran_ARCHIVE_FILE   netcdf-fortran-${NetCDF_Fortran_VERSION}.tar.gz)
set(NetCDF_Fortran_MD5_SUM        cc3bf530223e8f4aff93793b9f197bf3) 

#
# TPL: ExodusII
#
set(ExodusII_VERSION_MAJOR 6)
set(ExodusII_VERSION_MINOR 06)
set(ExodusII_VERSION ${ExodusII_VERSION_MAJOR}.${ExodusII_VERSION_MINOR})
set(ExodusII_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(ExodusII_ARCHIVE_FILE  exodus-${ExodusII_VERSION}.tar.gz)
set(ExodusII_MD5_SUM       cfd240dbc1251b08fb1d0ee2de40a44c)

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 2)
set(MSTK_VERSION_MINOR 23)
set(MSTK_VERSION_PATCH )
set(MSTK_VERSION ${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${MSTK_VERSION_PATCH})
set(MSTK_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MSTK_ARCHIVE_FILE   mstk-${MSTK_VERSION}.tgz)
set(MSTK_MD5_SUM        d1ceedfd43e18f8b5eee85b2d8e7d4fe)

#
# TPL: MOAB
#
set(MOAB_VERSION_MAJOR  r4276)
set(MOAB_VERSION_MINOR  )
set(MOAB_VERSION_PATCH  )
set(MOAB_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MOAB_ARCHIVE_FILE   MOAB-${MOAB_VERSION}.tar.gz)
set(MOAB_MD5_SUM        49da04e8905f6d730d92521e7ca7400e) 

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 11)
set(Trilinos_VERSION_MINOR 6)
set(Trilinos_VERSION_PATCH 1)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}.${Trilinos_VERSION_MINOR}.${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Trilinos_ARCHIVE_FILE   trilinos-${Trilinos_VERSION}-Source.tar.bz2)
set(Trilinos_MD5_SUM        b97d882535fd1856599b1c7338f5b45a)

#
# TPL: SEACAS
#  SEACAS is available in Trilinos 10.8 and above
set(SEACAS_VERSION_MAJOR 11)
set(SEACAS_VERSION_MINOR 6)
set(SEACAS_VERSION_PATCH 1)
set(SEACAS_VERSION ${SEACAS_VERSION_MAJOR}.${SEACAS_VERSION_MINOR}.${SEACAS_VERSION_PATCH})
set(SEACAS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SEACAS_ARCHIVE_FILE   trilinos-${SEACAS_VERSION}-Source.tar.bz2)
set(SEACAS_MD5_SUM        b97d882535fd1856599b1c7338f5b45a)

