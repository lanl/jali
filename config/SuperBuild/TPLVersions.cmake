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
#                   Xerces          3.1.1
#                   MPI version     1.4.4
#                   Curl            7.37.0
#                   Zlib            1.2.6
#                   Metis           5.1.0
#                   Boost           1.57
#                   BoostCmake      1.46.1
#                   HDF5            1.8.8
#                   NetCDF          4.3.2
#                   NetCDF fortran  4.2
#                   Exodus II       6.06
#                   MSTK            2.25
#                   Trilinos        12.2.2
#                   Seacas          12.2.2
#   1.0.1       - Boost updated to 1.58
#   1.0.2       - MSTK updated to 2.25rc1
#   1.0.3       - Trilinos updated to 12.2.2
#   1.0.4       - MSTK updated to 2.25
#   1.0.5       - MSTK updated to 2.26rc1 and Jali-tpl-config.cmake includes variables indicating if METIS and/or ZOLTAN are enabled
#   1.0.6       - MSTK updated to 2.26rc2
#   1.0.7       - MSTK updated to 2.27rc2
#   1.0.8       - MSTK updated to 2.27rc3
#   1.0.9       - MSTK updated to 2.27 and UnitTest++ updated to 1.6.0
#   1.1.0       - MSTK updated to 3.0.3, UnitTest++ to 2.0.0, Trilinos to 12.10.1, HDF5 1.8.18, Seacas to #173a1e6, zlib to 1.2.11, Boost to 1.63.0, NetCDF to 4.5.0, MOAB to 5.0.0, 

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

#
# Default location of TPLs on GitHub (poach Amanzi's TPL collection)
#
set (JALI_TPLS_DOWNLOAD_URL "https://raw.githubusercontent.com/amanzi/amanzi-tpls/master/src")

set (JALI_TPLS_VERSION_MAJOR 1)
set (JALI_TPLS_VERSION_MINOR 1)
set (JALI_TPLS_VERSION_PATCH 0)
set (JALI_TPLS_VERSION ${JALI_TPLS_VERSION}.${JALI_TPLS_VERSION_MINOR}.${JALI_TPLS_VERSION_PATCH})
#   Not sure how to create a meaningful hash key for the collection

#
# TPL: OpenMPI
#
set(OpenMPI_VERSION_MAJOR 2)
set(OpenMPI_VERSION_MINOR 1)
set(OpenMPI_VERSION_PATCH 2)
set(OpenMPI_VERSION ${OpenMPI_VERSION_MAJOR}.${OpenMPI_VERSION_MINOR}.${OpenMPI_VERSION_PATCH})
set(OpenMPI_URL_STRING     "https://www.open-mpi.org/software/ompi/v2.1/downloads/")
set(OpenMPI_ARCHIVE_FILE   openmpi-${OpenMPI_VERSION}.tar.bz2)
set(OpenMPI_SAVEAS_FILE    ${OpenMPI_ARCHIVE_FILE})
set(OpenMPI_MD5_SUM        ff2e55cc529802e7b0738cf87acd3ee4)
#
# TPL: CURL
#
set(CURL_VERSION_MAJOR 7)
set(CURL_VERSION_MINOR 56)
set(CURL_VERSION_PATCH 1)
set(CURL_VERSION ${CURL_VERSION_MAJOR}_${CURL_VERSION_MINOR}_${CURL_VERSION_PATCH})
set(CURL_URL_STRING     "https://github.com/curl/curl/archive")
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION_MAJOR}_${CURL_VERSION_MINOR}_${CURL_VERSION_PATCH}.tar.gz)
set(CURL_SAVEAS_FILE    curl-${CURL_VERSION}.tar.gz)
set(CURL_MD5_SUM        48c0db0d7b1407e19c51e8ef4f798d78)

#
# TPL: zlib
#
set(ZLIB_VERSION_MAJOR 1)
set(ZLIB_VERSION_MINOR 2)
set(ZLIB_VERSION_PATCH 11)
set(ZLIB_VERSION ${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH})
set(ZLIB_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(ZLIB_ARCHIVE_FILE   zlib-${ZLIB_VERSION}.tar.gz)
set(ZLIB_SAVEAS_FILE    ${ZLIB_ARCHIVE_FILE})
set(ZLIB_MD5_SUM        1c9f62f0778697a09d36121ead88e08e) 

#
# TPL: METIS
#
set(METIS_VERSION_MAJOR 5)
set(METIS_VERSION_MINOR 1)
set(METIS_VERSION_PATCH 0)
set(METIS_VERSION ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
set(METIS_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(METIS_ARCHIVE_FILE   metis-${METIS_VERSION}.tar.gz)
set(METIS_SAVEAS_FILE    ${METIS_ARCHIVE_FILE})
set(METIS_MD5_SUM        5465e67079419a69e0116de24fce58fe)

#
# TPL: UnitTest
#
set(UnitTest_VERSION_MAJOR 2)
set(UnitTest_VERSION_MINOR 0)
set(UnitTest_VERSION_PATCH 0)
set(UnitTest_VERSION ${UnitTest_VERSION_MAJOR}.${UnitTest_VERSION_MINOR}.${UnitTest_VERSION_PATCH})
set(UnitTest_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(UnitTest_ARCHIVE_FILE   unittest-cpp-${UnitTest_VERSION}.tgz)
set(UnitTest_SAVEAS_FILE    ${UnitTest_ARCHIVE_FILE})
set(UnitTest_MD5_SUM      29f958e355e516e7ab016b467974728d) 

#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 63)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.bz2)
set(Boost_SAVEAS_FILE    ${Boost_ARCHIVE_FILE})
set(Boost_MD5_SUM        1c837ecd990bb022d07e7aab32b09847)

#
# TPL: HDF5
#
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 18)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(HDF5_ARCHIVE_FILE   hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_SAVEAS_FILE    ${HDF5_ARCHIVE_FILE})
set(HDF5_MD5_SUM        dd2148b740713ca0295442ec683d7b1c)


#
# TPL: NetCDF
#
set(NetCDF_VERSION_MAJOR 4)
set(NetCDF_VERSION_MINOR 5)
set(NetCDF_VERSION_PATCH 0)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     "https://github.com/Unidata/netcdf-c/archive/")
set(NetCDF_ARCHIVE_FILE   v${NetCDF_VERSION}.tar.gz)
set(netCDF_SAVEAS_FILE    netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        a523ad253bd832efa632847940c2317e)

#
# TPL: NetCDF Fortran
#
set(NetCDF_Fortran_VERSION_MAJOR 4)
set(NetCDF_Fortran_VERSION_MINOR 4)
set(NetCDF_Fortran_VERSION_PATCH 4)
set(NetCDF_Fortran_VERSION ${NetCDF_Fortran_VERSION_MAJOR}.${NetCDF_Fortran_VERSION_MINOR}.${NetCDF_Fortran_VERSION_PATCH})
set(NetCDF_Fortran_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(NetCDF_Fortran_URL_STRING     "https://github.com/Unidata/netcdf-fortran/archive/")
set(NetCDF_Fortran_ARCHIVE_FILE   v${NetCDF_Fortran_VERSION}.tar.gz)
set(NetCDF_Fortran_SAVEAS_FILE    netcdf-fortran-${NetCDF_Fortran_VERSION}.tar.gz)
set(NetCDF_Fortran_MD5_SUM        418c7e998e63e6d76b2da14019fa9c8f) 

#
# TPL: ExodusII
#
set(ExodusII_VERSION_MAJOR 6)
set(ExodusII_VERSION_MINOR 06)
set(ExodusII_VERSION ${ExodusII_VERSION_MAJOR}.${ExodusII_VERSION_MINOR})
set(ExodusII_URL_STRING    ${JALI_TPLS_DOWNLOAD_URL})
set(ExodusII_ARCHIVE_FILE  exodus-${ExodusII_VERSION}.tar.gz)
set(ExodusII_MD5_SUM       cfd240dbc1251b08fb1d0ee2de40a44c)

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 3)
set(MSTK_VERSION_MINOR 0)
set(MSTK_VERSION_PATCH 3)
if (MSTK_VERSION_PATCH)
  set(MSTK_VERSION ${MSTK_VERSION_MAJOR}_${MSTK_VERSION_MINOR}_${MSTK_VERSION_PATCH})
else (MSTK_VERSION_PATCH)
  set(MSTK_VERSION ${MSTK_VERSION_MAJOR}_${MSTK_VERSION_MINOR})
endif (MSTK_VERSION_PATCH)
set(MSTK_URL_STRING     "https://github.com/MeshToolkit/MSTK/archive")
set(MSTK_ARCHIVE_FILE   v${MSTK_VERSION}.tar.gz)
set(MSTK_MD5_SUM        a9561f1175da2c2863c517b08abe5f16)

#
# TPL: MOAB
#

set(MOAB_VERSION_MAJOR  5)
set(MOAB_VERSION_MINOR  0)
set(MOAB_VERSION_PATCH  0)
set(MOAB_URL_STRING     "ftp://ftp.mcs.anl.gov/pub/fathom")
set(MOAB_ARCHIVE_FILE   MOAB-${MOAB_VERSION}.tar.gz)
set(MOAB_MD5_SUM        1840ca02366f4d3237d44af63e239e3b)

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 12)
set(Trilinos_VERSION_MINOR 10)
set(Trilinos_VERSION_PATCH 1)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}-${Trilinos_VERSION_MINOR}-${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING      "https://github.com/trilinos/Trilinos/archive")
set(Trilinos_ARCHIVE_FILE   trilinos-release-${Trilinos_VERSION}.tar.gz)
set(Trilinos_MD5_SUM        667333dbd7c0f031d47d7c5511fd0810)

#
# TPL: SEACAS
set(SEACAS_VERSION_MAJOR 173a1e6)
set(SEACAS_VERSION_MINOR 0)
set(SEACAS_VERSION_PATCH 0)
set(SEACAS_VERSION ${SEACAS_VERSION_MAJOR})
set(SEACAS_URL_STRING     ${JALI_TPLS_DOWNLOAD_URL})
set(SEACAS_ARCHIVE_FILE   seacas-${SEACAS_VERSION}.tgz)
set(SEACAS_SAVEAS_FILE    ${SEACAS_ARCHIVE_FILE})
set(SEACAS_MD5_SUM        3235d1b885ee8e1a04408382f50bd0f0)

