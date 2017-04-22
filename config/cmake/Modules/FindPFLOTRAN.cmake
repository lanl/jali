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
# Jali PFLOTRAN Find Module
#
# Usage:
#    Control the search through PFLOTRAN_DIR or setting environment variable
#    PFLOTRAN_ROOT to the PFLOTRAN installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PFLOTRAN_FOUND            (BOOL)       Flag indicating if PFLOTRAN was found
#    PFLOTRAN_INCLUDE_DIR      (PATH)       Path to the PFLOTRAN include file
#    PFLOTRAN_INCLUDE_DIRS     (LIST)       List of all required include files
#    PFLOTRAN_LIBRARY_DIR      (PATH)       Path to the PFLOTRAN library
#    PFLOTRAN_LIBRARY          (FILE)       PFLOTRAN library
#    PFLOTRAN_LIBRARIES        (LIST)       List of all required PFLOTRAN libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Jali CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( PFLOTRAN_LIBRARIES AND PFLOTRAN_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(PFLOTRAN_LIBRARIES AND PFLOTRAN_INCLUDE_DIRS)

    # Cache variables
    if(PFLOTRAN_DIR)
      set(PFLOTRAN_DIR "${PFLOTRAN_DIR}" CACHE PATH "Path to search for PFLOTRAN include and library files")
    endif()

    if(PFLOTRAN_INCLUDE_DIR)
      set(PFLOTRAN_INCLUDE_DIR "${PFLOTRAN_INCLUDE_DIR}" CACHE PATH "Path to search for PFLOTRAN include files")
    else()
      find_path(PFLOTRAN_INCLUDE_DIR pflotran_alquimia_interface.h ${PFLOTRAN_DIR}/include)
      if ( NOT PFLOTRAN_INCLUDE_DIR )
        message(SEND_ERROR "Cannot locate PFLOTRAN include directory")
      endif()
    endif()

    if(PFLOTRAN_LIBRARY_DIR)
        set(PFLOTRAN_LIBRARY_DIR "${PFLOTRAN_LIBRARY_DIR}" CACHE PATH "Path to search for PFLOTRAN library files")
    else()
      find_path(PFLOTRAN_LIBRARY_DIR NAMES libpflotranchem.a PATHS ${PFLOTRAN_DIR}/lib)
      if ( NOT PFLOTRAN_LIBRARY_DIR )
        message(SEND_ERROR "Cannot locate PFLTORAN library directory")
      endif()
    endif()

    # Search for libraries 
    
    set( PFLOTRAN_TARGET pflotranchem )

    find_library(_PFLOTRAN_LIBRARY
                 NAMES pflotranchem
                 PATHS ${PFLOTRAN_LIBRARY_DIR}
                 NO_DEFAULT_PATH )

    if ( _PFLOTRAN_LIBRARY ) 
      add_imported_library(${PFLOTRAN_TARGET}
	                   LOCATION ${_PFLOTRAN_LIBRARY}
                           LINK_LANGUAGES "Fortran")
      set(PFLOTRAN_LIBRARY ${PFLOTRAN_TARGET})
    else()
      set(PFLOTRAN_LIBRARY PFLOTRAN_LIBRARY-NOTFOUND)
      message(SEND_ERROR "Cannot locate PFLOTRAN library")
    endif()    
    
   
    # Define the LIBRARIES and INCLUDE_DIRS

    set(PFLOTRAN_INCLUDE_DIRS ${PFLOTRAN_INCLUDE_DIR})
    set(PFLOTRAN_LIBRARIES    ${PFLOTRAN_LIBRARY})

endif(PFLOTRAN_LIBRARIES AND PFLOTRAN_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(PFLOTRAN DEFAULT_MSG
                                           PFLOTRAN_INCLUDE_DIRS
					   PFLOTRAN_LIBRARIES)

# find_package)handle)standard_args should set PFLOTRAN_FOUND but it does not!
if ( PFLOTRAN_LIBRARIES AND PFLOTRAN_INCLUDE_DIRS)
    set(PFLOTRAN_FOUND TRUE)
else()
    set(PFLOTRAN_FOUND FALSE)
endif()

mark_as_advanced(
  PFLOTRAN_INCLUDE_DIR
  PFLOTRAN_INCLUDE_DIRS
  PFLOTRAN_LIBRARY
  PFLOTRAN_LIBRARIES
  PFLOTRAN_LIBRARY_DIR
)
