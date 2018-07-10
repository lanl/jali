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
# CHECK_DOWNLOAD
#
# USAGE:
#  CHECK_DOWNLOAD(TEST_URL web-address
#                 TEST_FILE file
#                 TIMEOUT seconds)
#
# Attempt to download file from web-address using the system curl binary.
# Default TIMEOUT for this command is 60 seconds. If download fails, throws
# fatal error. 
#
include(CMakeParseArguments)
include(PrintVariable)

function(CHECK_DOWNLOAD)


  # Parse the arguments
  set(_flags    "")
  set(_oneValue "TEST_URL;TEST_FILE;TIMEOUT")
  set(_multiValue "")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})

  # Default timeout is 60 seconds
  set(command_timeout 60)
  if ( PARSE_TIMEOUT )
    set(command_timeout ${PARSE_TIMEOUT})
  endif()  

  # Need a URL AND FILE name
  if ( PARSE_TEST_URL AND PARSE_TEST_FILE )
    set(url_string ${PARSE_TEST_URL}/${PARSE_TEST_FILE})
  else()
    message(FATAL_ERROR "Invalid arguments to CHECK_DOWNLOAD. "
                        "Must define TEST_URL AND TEST_FILE")
  endif()                    
       
  message(STATUS "Checking external downloads (${url_string})")

  file(DOWNLOAD
       ${url_string}
       ${CMAKE_CURRENT_BINARY_DIR}/${PARSE_TEST_FILE}
       TIMEOUT ${command_timeout}
       STATUS result)
  list(GET result 0 ret_code)
  list(GET result 1 error_str)
  if ( "${ret_code}" EQUAL 0 )
    message(STATUS "Checking external downloads (${url_string}) -- works")
  else()  
    message(SEND_ERROR "Failed to download ${url_string} "
                       "Return Code:${ret_code}\n"
                       "Output:\n${error_str}\n"
                       "You can disable external downloads with"
                       "\n-DDISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE\n"
                       "If external downloads are disabled, ALL TPL source files must be found in one directory"
                       " and define this directory with"
                       "\n-D TPL_DOWNLOAD_DIR:FILEPATH=\n")
    message(FATAL_ERROR "Failed to download ${PARSE_TEST_FILE} from ${PARSE_TEST_URL}")

  endif()  


endfunction(CHECK_DOWNLOAD)
