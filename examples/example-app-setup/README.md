# How to set up an app to use Jali

This directory has a simple example of how to setup a CMake build system to compile your app with Jali

The source for your app is assumed to be under 'src'
Jali_TOP_DIR is under which Jali lib and include dirs are installed

So here are the steps:

   module load some_compiler (e.g. intel/15.0.3)
   module load some_openmpi_wrappers (e.g. openmpi/1.6.5)
   SOURCE=`pwd`
   JALI_TOP_DIR=/Path/to/where/Jali/includes/and/libs/are/installed
   mkdir build
   cd build
   cmake \
       -D CMAKE_C_COMPILER:FILEPATH=`which mpicc` \
       -D CMAKE_CXX_COMPILER:FILEPATH=`which mpiCC` \
       -D Jali_DIR:FILEPATH=${Jali_TOP_DIR}/lib \
       ${SOURCE}
   make

The top level CMakeLists.txt file is setup to pull in the Jali includes and libraries from the variable Jali_DIR defined in the cmake command

You can get more sophisticated with the CMAKE setup as your project develops or your configuration changes. For example you can put the cmake command in a do-configure script and use that to configure the project in the build directory before running 'make'. 

