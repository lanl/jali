# Installation instructions

There are a few configuration examples in `config/configure-examples/`. Below
we list copy & paste instructions for several machines. You can easily adapt
them for other machines.

## Darwin

Execute the following from the Jali root directory:

    module load compilers/gcc/4.9.2 mpi/openmpi-1.8.4-intel_15.0.3 cmake
    SOURCE=`pwd`
    TPL_INSTALL_PREFIX=$SOURCE/inst-tpl
    JALI_INSTALL_PREFIX=$SOURCE/inst-jali
    mkdir build-tpl
    cd build-tpl
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_Fortran_COMPILER=`which mpif90` \
        -D DISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE \
        -D TPL_DOWNLOAD_DIR:PATH=/usr/projects/ngc/public/tpl-downloads/ \
        -D TPL_INSTALL_PREFIX=$TPL_INSTALL_PREFIX \
        $SOURCE/config/SuperBuild/
    make -j16
    make install
    cd ..
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Debug \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j16
    ctest -j16
    make install

## Varan

Execute the following from the Jali root directory:

    export MODULEPATH=""
    . /opt/local/packages/Modules/default/init/sh
    module load intel/15.0.3 openmpi/1.6.5
    SOURCE=`pwd`
    TPL_INSTALL_PREFIX=$SOURCE/inst-tpl
    JALI_INSTALL_PREFIX=$SOURCE/inst-jali
    mkdir build-tpl
    cd build-tpl
    cmake \
        -D CMAKE_C_COMPILER=`which mpicc` \
        -D CMAKE_CXX_COMPILER=`which mpiCC` \
        -D CMAKE_Fortran_COMPILER=`which mpif90` \
        -D DISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE \
        -D TPL_DOWNLOAD_DIR:PATH=/usr/local/codes/ngc/private/tpl-downloads2/ \
        -D TPL_INSTALL_PREFIX=$TPL_INSTALL_PREFIX \
        $SOURCE/config/SuperBuild/
    make -j16
    make install
    cd ..
    mkdir build-jali
    cd build-jali
    cmake \
      -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake \
      -D CMAKE_BUILD_TYPE=Debug \
      -D CMAKE_CXX_FLAGS='-std=c++11' \
      -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX \
      -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE \
      -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
      -D ENABLE_MSTK_Mesh:BOOL=TRUE \
      -D ENABLE_STK_Mesh:BOOL=FALSE \
      -D ENABLE_MOAB_Mesh:BOOL=FALSE \
      ${SOURCE}
    make -j16
    ctest -j16
    make install
