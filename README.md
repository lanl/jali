#
#    Jali VERSION 0.6.5 (Copyright, Los Alamos National Laboratory)
#

# Summary

Jali is a parallel unstructured mesh infrastructure library designed
for use by multi-physics simulations. It supports 2D and 3D arbitrary
polyhedral meshes distributed over hundreds to thousands of
nodes. Jali can read and write Exodus II meshes along with fields and sets
on the mesh. Jali is built upon MSTK, an open source general purpose
unstructured mesh infrastructure library from LANL. Jali is
copyrighted as open source but is not yet available to the public. 

# Third Party Libraries

Jali uses a number of third party libraries (TPLs) primarily for mesh
import/export, mesh partitioning and MPI communication. To build Jali,
you must point the build system to a directory containing the TPLs or
build the third party libraries as described below.

Jali *0.9.2*  uses version *1.0.7* or higher of the TPL set. See
$JALI_SOURCE/config/SuperBuild/TPLVersions.cmake for details.

# Installation instructions

There are a few configuration examples in
`config/configure-examples/`. Below we list copy & paste instructions
for several machines. The first option is the *easy* option in which
pre-installed third-party libraries (TPLs) are used to build Jali. The
second option is for cases where you want to build all the TPLs used
by Jali and build Jali using these custom TPLs. Naturally, in the
following instructions, you may use your choice of compiler, openmpi
libraries and cmake version, subject to compatibility requirements of Jali


## Darwin

If you want to build only Jali, execute the following from the Jali
root directory (you may have to modify the exact version number of the
TPLs directory depending on how out-of-sync this file is with whats on
your system):

    /bin/tcsh
    module load openmpi/1.10.0-intel_15.0.3 cmake
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/projects/ngc/private/jali-1.0.4-tpl-intel
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit

If you want to build a custom set of Jali TPLs and build Jali using
these TPLs, execute the following from the Jali root directory: 

    /bin/tcsh
    module load openmpi/1.10.0-intel_15.0.3 cmake
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX $SOURCE/inst-tpl
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-tpl
    cd build-tpl
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_Fortran_COMPILER=`which mpif90` \
        -D DISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE \
        -D TPL_DOWNLOAD_DIR:PATH=/usr/projects/ngc/private/tpl-downloads/ \
        -D TPL_INSTALL_PREFIX=$TPL_INSTALL_PREFIX \
        $SOURCE/config/SuperBuild/
    make -j
    make install
    cd ..
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit

## Moonlight (or Pinto)

If you want to build only Jali, execute the following from the Jali
root directory (you may have to modify the exact version number of the
TPLs directory depending on how out-of-sync this file is with whats on
your system):

    /bin/tcsh
    module load intel/15.0.3 openmpi/1.6.5
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/projects/ngc/private/jali-1.0.4-tpl-intel
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit

If you want to build a custom set of Jali TPLs and build Jali using
these TPLs, execute the following from the Jali root directory: 

    /bin/tcsh
    module load intel/15.0.3 openmpi/1.6.5
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX $SOURCE/inst-tpl
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-tpl
    cd build-tpl
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_Fortran_COMPILER=`which mpif90` \
        -D DISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE \
        -D TPL_DOWNLOAD_DIR:PATH=/usr/projects/ngc/private/tpl-downloads/ \
        -D TPL_INSTALL_PREFIX=$TPL_INSTALL_PREFIX \
        $SOURCE/config/SuperBuild/
    make -j
    make install
    cd ..
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit

## XLAN (Varan or Barugon)

If you want to build only Jali, execute the following from the Jali
root directory (you may have to modify the exact version number of the
TPLs directory depending on how out-of-sync this file is with whats on
your system):

    /bin/tcsh
    /opt/local/packages/Modules/default/init/sh
    module load intel/15.0.3 openmpi/1.6.5
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/local/codes/ngc/private/jali-1.0.4-tpl-intel
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit

If you want to build a custom set of Jali TPLs and build Jali using
these TPLs, execute the following from the Jali root directory: 

    /bin/tcsh
    /opt/local/packages/Modules/default/init/sh
    module load intel/15.0.3 openmpi/1.6.5
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX $SOURCE/inst-tpl
    setenv JALI_INSTALL_PREFIX $SOURCE/inst-jali
    mkdir build-tpl
    cd build-tpl
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_Fortran_COMPILER=`which mpif90` \
        -D DISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE \
        -D TPL_DOWNLOAD_DIR:PATH=/usr/local/codes/ngc/private/tpl-downloads/ \
        -D TPL_INSTALL_PREFIX=$TPL_INSTALL_PREFIX \
        $SOURCE/config/SuperBuild/
    make -j
    make install
    cd ..
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j
    make test
    make install
    exit
