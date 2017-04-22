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
# Jali ASCEMIO Find Module
#
# Usage:
#    Control the search through ASCEMIO_DIR or setting environment variable
#    ASCEMIO_ROOT to the ASCEMIO installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    ASCEMIO_FOUND            (BOOL)       Flag indicating if ASCEMIO was found
#    ASCEMIO_INCLUDE_DIR      (PATH)       Path to the ASCEMIO include file
#    ASCEMIO_INCLUDE_DIRS     (LIST)       List of all required include files
#    ASCEMIO_LIBRARY_DIR      (PATH)       Path to the ASCEMIO library
#    ASCEMIO_LIBRARY          (FILE)       ASCEMIO library
#    ASCEMIO_LIBRARIES        (LIST)       List of all required ASCEMIO libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Jali CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( ASCEMIO_LIBRARIES AND ASCEMIO_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(ASCEMIO_LIBRARIES AND ASCEMIO_INCLUDE_DIRS)

    # Cache variables
    if(ASCEMIO_DIR)
        set(ASCEMIO_DIR "${ASCEMIO_DIR}" CACHE PATH "Path to search for ASCEMIO include and library files")
    endif()

    if(ASCEMIO_INCLUDE_DIR)
        set(ASCEMIO_INCLUDE_DIR "${ASCEMIO_INCLUDE_DIR}" CACHE PATH "Path to search for ASCEMIO include files")
    endif()

    if(ASCEMIO_LIBRARY_DIR)
        set(ASCEMIO_LIBRARY_DIR "${ASCEMIO_LIBRARY_DIR}" CACHE PATH "Path to search for ASCEMIO library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) ASCEMIO_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) ASCEMIO_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(ascemio_inc_names "parallelIO.h")
    if (ASCEMIO_INCLUDE_DIR)

        if (EXISTS "${ASCEMIO_INCLUDE_DIR}")

            find_path(io_test_include_path
                      NAMES ${ascemio_inc_names}
                      HINTS ${ASCEMIO_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT io_test_include_path)
                message(SEND_ERROR "Can not locate ${netcdf_inc_names} in ${ASCEMIO_INCLUDE_DIR}")
            endif()
            set(ASCEMIO_INCLUDE_DIR "${io_test_include_path}")

        else()
            message(SEND_ERROR "ASCEMIO_INCLUDE_DIR=${ASCEMIO_INCLUDE_DIR} does not exist")
            set(ASCEMIO_INCLUDE_DIR "ASCEMIO_INCLUDE_DIR-NOTFOUND")
        endif()

    else() 

	set(ascemio_inc_suffixes "include")
        if(ASCEMIO_DIR)

            if (EXISTS "${ASCEMIO_DIR}" )

		find_path(ASCEMIO_INCLUDE_DIR
		          NAMES ${ascemio_inc_names}
		          HINTS ${ASCEMIO_DIR}
		          PATH_SUFFIXES ${ascemio_inc_suffixes}
		          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ASCEMIO_DIR=${ASCEMIO_DIR} does not exist")
                 set(ASCEMIO_INCLUDE_DIR "ASCEMIO_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(ASCEMIO_INCLUDE_DIR
                      NAMES ${ascemio_inc_names}
                      PATH_SUFFIXES ${ascemio_inc_suffixes})

        endif()

    endif()


    if ( NOT ASCEMIO_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate ASCEMIO include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) ASCEMIO_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) ASCEMIO_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    if (ASCEMIO_LIBRARY_DIR)

        if (EXISTS "${ASCEMIO_LIBRARY_DIR}")

            find_library(_ASCEMIO_LIBRARY
                         NAMES parallelio
                         HINTS ${ASCEMIO_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "ASCEMIO_LIBRARY_DIR=${ASCEMIO_LIBRARY_DIR} does not exist")
            set(ASCEMIO_LIBRARY "ASCEMIO_LIBRARY-NOTFOUND")
        endif()

    else() 

        if(ASCEMIO_DIR)

            if (EXISTS "${ASCEMIO_DIR}" )

                find_library(_ASCEMIO_LIBRARY
                             NAMES parallelio
                             HINTS ${ASCEMIO_DIR}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ASCEMIO_DIR=${ASCEMIO_DIR} does not exist")
                 set(ASCEMIO_LIBRARY "ASCEMIO_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(_ASCEMIO_LIBRARY
                         NAMES parallelio
                         PATH_SUFFIXES ${ascemio_lib_suffixes})

        endif()

    endif()

    if ( _ASCEMIO_LIBRARY )
        set(ASCEMIO_LIBRARY ascemio)
        add_imported_library(${ASCEMIO_LIBRARY}
	                     LOCATION ${_ASCEMIO_LIBRARY}
                             LINK_LANGUAGES "C")
    else()
        set(ASCEMIO_LIBRARY ASCEMIO_LIBRARY-NOTFOUND)
        message(SEND_ERROR "Can not locate ASCEMIO library")
    endif()    
    
   
    # Define the LIBRARIES and INCLUDE_DORS
    set(ASCEMIO_INCLUDE_DIRS ${ASCEMIO_INCLUDE_DIR})
    set(ASCEMIO_LIBRARIES    ${ASCEMIO_LIBRARY})

    message(STATUS "ASCEMIO requires HDF5")
    find_package(HDF5 QUIET REQUIRED)
    set_target_properties(${ASCEMIO_LIBRARY} PROPERTIES
                          IMPORTED_LINK_INTERFACE_LIBRARIES "${HDF5_C_LIBRARIES}")
    list(APPEND ASCEMIO_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})			
    #add_package_dependency(ASCEMIO DEPENDS_ON HDF5)

endif(ASCEMIO_LIBRARIES AND ASCEMIO_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(ASCEMIO DEFAULT_MSG
                                           ASCEMIO_INCLUDE_DIRS
					   ASCEMIO_LIBRARIES)

# find_package)handle)standard_args should set ASCEMIO_FOUND but it does not!
if ( ASCEMIO_LIBRARIES AND ASCEMIO_INCLUDE_DIRS)
    set(ASCEMIO_FOUND TRUE)
else()
    set(ASCEMIO_FOUND FALSE)
endif()

mark_as_advanced(
  ASCEMIO_INCLUDE_DIR
  ASCEMIO_INCLUDE_DIRS
  ASCEMIO_LIBRARY
  ASCEMIO_LIBRARIES
  ASCEMIO_LIBRARY_DIR
)
