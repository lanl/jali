# Copyright (c) 2017, Los Alamos National Security, LLC
# All rights reserved.

# Copyright 2017. Los Alamos National Security, LLC. This software was
# produced under U.S. Government contract DE-AC52-06NA25396 for Los
# Alamos National Laboratory (LANL), which is operated by Los Alamos
# National Security, LLC for the U.S. Department of Energy. The
# U.S. Government has rights to use, reproduce, and distribute this
# software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
# LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
# derivative works, such modified software should be clearly marked, so
# as not to confuse it with the version available from LANL.
 
# Additionally, redistribution and use in source and binary forms, with
# or without modification, are permitted provided that the following
# conditions are met:

# 1.  Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2.  Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 3.  Neither the name of Los Alamos National Security, LLC, Los Alamos
# National Laboratory, LANL, the U.S. Government, nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
 
# THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
# ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#
# Jali MOAB Find Module
#
# Usage:
#    Control the search through MOAB_DIR or setting environment variable
#    MOAB_ROOT to the MOAB installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    MOAB_FOUND            (BOOL)       Flag indicating if MOAB was found
#    MOAB_INCLUDE_DIR      (PATH)       Path to the MOAB include file
#    MOAB_INCLUDE_DIRS     (LIST)       List of all required include files
#    MOAB_LIBRARY_DIR      (PATH)       Path to the MOAB library
#    MOAB_LIBRARY          (FILE)       MOAB library
#    MOAB_LIBRARIES        (LIST)       List of all required MOAB libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Jali CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)
include(AddImportedLibrary)

if ( MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)

    # Cache variables
    if(MOAB_DIR)
        set(MOAB_DIR "${MOAB_DIR}" CACHE PATH "Path to search for MOAB include and library files")
    endif()

    if(MOAB_INCLUDE_DIR)
        set(MOAB_INCLUDE_DIR "${MOAB_INCLUDE_DIR}" CACHE PATH "Path to search for MOAB include files")
    endif()

    if(MOAB_LIBRARY_DIR)
        set(MOAB_LIBRARY_DIR "${MOAB_LIBRARY_DIR}" CACHE PATH "Path to search for MOAB library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) MOAB_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) MOAB_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(moab_inc_names "MBCore.hpp")
    if (MOAB_INCLUDE_DIR)

        if (EXISTS "${MOAB_INCLUDE_DIR}")

            find_path(moab_test_include_path
                      NAMES ${moab_inc_names}
                      HINTS ${MOAB_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT moab_test_include_path)
                message(SEND_ERROR "Can not locate ${moab_inc_names} in ${MOAB_INCLUDE_DIR}")
            endif()
            set(MOAB_INCLUDE_DIR "${moab_test_include_path}")

        else()
            message(SEND_ERROR "MOAB_INCLUDE_DIR=${MOAB_INCLUDE_DIR} does not exist")
            set(MOAB_INCLUDE_DIR "MOAB_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(moab_inc_suffixes "include")
        if(MOAB_DIR)

            if (EXISTS "${MOAB_DIR}" )

                find_path(MOAB_INCLUDE_DIR
                          NAMES ${moab_inc_names}
                          HINTS ${MOAB_DIR}
                          PATH_SUFFIXES ${moab_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "MOAB_DIR=${MOAB_DIR} does not exist")
                 set(MOAB_INCLUDE_DIR "MOAB_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(MOAB_INCLUDE_DIR
                      NAMES ${moab_inc_names}
                      PATH_SUFFIXES ${moab_inc_suffixes})

        endif()

    endif()

    if ( NOT MOAB_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate MOAB include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) MOAB_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) MOAB_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(moab_lib_names "MOAB")
    if (MOAB_LIBRARY_DIR)

        if (EXISTS "${MOAB_LIBRARY_DIR}")

            find_library(_MOAB_LIBRARY
                         NAMES ${moab_lib_names}
                         HINTS ${MOAB_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "MOAB_LIBRARY_DIR=${MOAB_LIBRARY_DIR} does not exist")
            set(_MOAB_LIBRARY "_MOAB_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND moab_lib_suffixes "lib" "Lib")
        if(MOAB_DIR)

            if (EXISTS "${MOAB_DIR}" )

                find_library(_MOAB_LIBRARY
                             NAMES ${moab_lib_names}
                             HINTS ${MOAB_DIR}
                             PATH_SUFFIXES ${moab_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "MOAB_DIR=${MOAB_DIR} does not exist")
                 set(_MOAB_LIBRARY "_MOAB_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(_MOAB_LIBRARY
                         NAMES ${moab_lib_names}
                         PATH_SUFFIXES ${moab_lib_suffixes})

        endif()

    endif()

    # Create the target
    if ( _MOAB_LIBRARY )
      set(MOAB_LIBRARY MOAB)
      add_imported_library(${MOAB_LIBRARY}
	                   LOCATION ${_MOAB_LIBRARY}
			   LINK_LANGUAGES "C;CXX")
    else()  			 
      message(SEND_ERROR "Can not locate MOAB library")
    endif()    

   
    # Update the INCLUDE_DIRS and LIBRARIES variables
    set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIR})
    set(MOAB_LIBRARIES    ${MOAB_LIBRARY})

    # Define the dependent libs
    set(_MOAB_DEP_LIBS)

    # HDF5
    find_package(HDF5 QUIET REQUIRED)
    list(APPEND MOAB_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    list(APPEND _MOAB_DEP_LIBS ${HDF5_LIBRARIES})

    # NetCDF (This should come with the HDF5)
    find_package(NetCDF QUIET REQUIRED)
    list(APPEND MOAB_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS})
    list(APPEND _MOAB_DEP_LIBS ${NetCDF_C_LIBRARIES})

    set_target_properties(${MOAB_LIBRARY} PROPERTIES
                          IMPORTED_LINK_INTERFACE_LIBRARIES "${_MOAB_DEP_LIBS}")

endif(MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(MOAB DEFAULT_MSG
                                           MOAB_LIBRARIES
                                           MOAB_INCLUDE_DIRS)

# Define the version

mark_as_advanced(
  MOAB_INCLUDE_DIR
  MOAB_INCLUDE_DIRS
  MOAB_LIBRARY
  MOAB_LIBRARIES
  MOAB_LIBRARY_DIR
)
