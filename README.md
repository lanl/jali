Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
National Laboratory, LANL, the U.S. Government, nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Summary

Jali is a parallel unstructured mesh infrastructure library designed
for use by multi-physics simulations. It supports 2D and 3D arbitrary
polyhedral meshes distributed over hundreds to thousands of
nodes. Jali can read and write Exodus II meshes along with fields and
sets on the mesh and support for other formats is partially
implemented or is in the plans. Jali is built upon MSTK
(https://github.com/MeshToolkit/MSTK), an open source general purpose
unstructured mesh infrastructure library from Los Alamos National
Laboratory. While it has been made to work with other mesh frameworks
such as MOAB and STKmesh in the past, support for maintaining the
interface to these frameworks has been suspended for now. Jali
supports distributed as well as on-node parallelism. Support of
on-node parallelism is through direct use of the mesh calls in
multi-threaded constructs or through use of "tiles" which are
submeshes or sub-partitions of a partition destined for a compute
node.

Jali is derived from the mesh infrastructure of the early versions of
the open source software Amanzi (whose three-part BSD open source
copyright assertion included Los Alamos National Laboratory, Pacific
Northwest National Laboratory and Lawrence Berkeley
Laboratory). Jali's copyright is being asserted solely under Los
Alamos National Laboratory having been rewritten to include only code
written by Los Alamos National Laboratory developers.

# Authors

Jali is derived from the unstructured mesh infrastructure in Amanzi (https://github.com/amanzi/amanzi) with non-LANL contributions stripped out. While much of the unstructured mesh infrastructure in  Amanzi (and subsequently Jali) was developed by Rao Garimella, its design was influenced by the entire Amanzi team. Additionally, contributions to the Amanzi or Jali code base have been made by Markus Berndt, Konstantin Lipnikov, Ethan Coon, Daniil Svyatsky, Chris Malone, Chris Sewell, Charles Ferenbaugh and Ondrej Certik. In addition, most of the build system in Jali was taken as-is from Amanzi. This build system was largely developed by Lori Pritchett-Sheats.

# Third Party Libraries

Jali uses a number of third party libraries (TPLs) primarily for mesh
import/export, mesh partitioning and MPI communication. To build Jali,
you must point the build system to a directory containing the TPLs or
build the third party libraries as described below.

Jali *1.0.0*  uses version *1.1.0* or higher of the TPL set. See
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
    module load openmpi/2.1.2-intel_17.0.6 cmake
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/projects/ngc/private/jali-tpl/1.1.0-intel-17.0.6-openmpi-2.1.2
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
    module load openmpi/2.1.2-intel_17.0.6 cmake
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

## Snow

If you want to build only Jali, execute the following from the Jali
root directory (you may have to modify the exact version number of the
TPLs directory depending on how out-of-sync this file is with whats on
your system):

    /bin/tcsh
    module load intel/17.0.4 openmpi/2.1.2
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/projects/ngc/private/jali-tpl/1.1.0-intel-17.0.4-openmpi-2.1.2
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
    module load intel/17.0.4 openmpi/2.1.2
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
    module load intel/17.0.1 openmpi/1.10.5
    setenv SOURCE `pwd`
    setenv TPL_INSTALL_PREFIX /usr/local/codes/ngc/private/jali-tpl/1.1.0-intel-17.0.1-openmpi-1.10.5
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
    module load intel/17.0.1 openmpi/1.10.5
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
