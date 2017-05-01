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
# Jali Mercurial Find Module
#
# Usage:
#
#   Control the search through the variable MERCURIAL_ROOT. Searches for the mercurial
#   binary (i.e. hg) and if it finds the binary defines the version
#
#
#  Module sets
#
#  MERCURIAL_EXECUTABLE   Full path name of the Mercurial binary
#  MERCURIAL_VERSION      Mercurial Version
#  MERCURIAL_FOUND        Flag indicating if mercurial has been found
#
#
include(FindPackageHandleStandardArgs)

if (MERCURIAL_NO_SYSTEM_PATHS)
  set(_hg_FIND_OPTIONS NO_CMAKE_SYSTEM_PATH)
endif()

if ( NOT MERCURIAL_EXECUTABLE )

  find_program(MERCURIAL_EXECUTABLE
               name hg
               PATHS ${MERCURIAL_ROOT}
               ${_hg_FIND_OPTIONS})

endif()

if ( MERCURIAL_EXECUTABLE )

  execute_process(COMMAND ${MERCURIAL_EXECUTABLE} --version
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE MERCURIAL_VERSION
                  ERROR_VARIABLE  err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

  if (exit_code)
    message(WARNING "Could not determine mercurial version:${err}")
  endif()

endif()

find_package_handle_standard_args(Mercurial DEFAULT_MSG 
                                  MERCURIAL_EXECUTABLE)
                                           


