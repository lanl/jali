# Copyright (c) 2019, Triad National Security, LLC
# All rights reserved.

# Copyright 2019. Triad National Security, LLC. This software was
# produced under U.S. Government contract 89233218CNA000001 for Los
# Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S. Department of Energy. 
# All rights in the program are reserved by Triad National Security,
# LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting
# on its behalf a nonexclusive, paid-up, irrevocable worldwide license
# in this material to reproduce, prepare derivative works, distribute
# copies to the public, perform publicly and display publicly, and to
# permit others to do so
 
# 
# This is open source software distributed under the 3-clause BSD license.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of Triad National Security, LLC, Los Alamos
#    National Laboratory, LANL, the U.S. Government, nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
 
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# ############################################################################ #
#
# DefineCompilerVersion
#  
# ############################################################################ #

include(PrintVariable)

function(DEFINE_COMPILER_VERSION)

   # Define languages
   set(_languages C CXX Fortran)
   foreach ( _lang IN LISTS _languages)

       # Define the full path name of the compiler and the id
       # Default value for each language is NOTFOUND
       set(_lang_compiler ${CMAKE_${_lang}_COMPILER})
       set(_compiler_id   ${CMAKE_${_lang}_COMPILER_ID})
       set(_version       CMAKE_${_lang}_COMPILER_VERSION-NOTFOUND)

       message(STATUS "Identifying ${_lang_compiler} (${_compiler_id}) version")

       # Check only if the ID is set
       if ( _compiler_id )
	     string(TOUPPER ${_compiler_id} _compiler_id_uc)
         set(_version_cmd_opt)
         set(_regexp_pattern)

         # For now, I assume the version option is the same for all languages in
         # a compiler group. This may need tweaking in the future.
         if ( ${_compiler_id_uc} STREQUAL "GNU" )
           set(_version_cmd_opt "--version")
           set(_regexp_pattern ".*\(GCC\)[ ]+([0-9]+\\.[0-9]+\\.[0-9]+).*")
         elseif(${_compiler_id_uc} STREQUAL "PGI" )
           set(_version_cmd_opt "-V")
           set(_regexp_pattern ".*pgcc[ ]+([0-9]+\\.[0-9]+-[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "INTEL" )
           set(_version_cmd_opt "-V")
	   set(_regexp_pattern ".*,[ ]+Version ([0-9]+\\.[0-9]+\\.*[0-9]*).*Build.*")
         elseif (${_compiler_id_uc} STREQUAL "PATHSCALE" )
           set(_version_cmd_opt "--version")
           set(_regexp_pattern ".*Version ([0-9]+\\.[0-9]+\\.[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "CRAY" )
           set(_version_cmd_opt "-V")
           set(_regexp_pattern ".*Cray[ ]+C[ ]+:[ ]+Version ([0-9]+\\.[0-9]+\\.[0-9]+).*")
         else()
           message(WARNING "Unknown compiler ID type ${_lang_id_var}")
         endif()

         # Execute the command if the option was set
         if (DEFINED _version_cmd_opt)
           execute_process(COMMAND ${_lang_compiler} ${_version_cmd_opt}
                           RESULT_VARIABLE result
                           OUTPUT_VARIABLE output
                           ERROR_VARIABLE  output 
                           OUTPUT_STRIP_TRAILING_WHITESPACE
                           ERROR_STRIP_TRAILING_WHITESPACE)
  
           # result > 0 indicates an error return 
           if ( result ) 
             message(SEND_ERROR "${_lang_compiler} ${_version_cmd_opt} failed."
                                "Ouptut:\n${output}")
           else()
	     string(REGEX REPLACE "${_regexp_pattern}" "\\1" _version "${output}")
             # A fix for the PGI case because they use a '-' to separate the patch number. Annoying.
             if ( ${_compiler_id_uc} STREQUAL "PGI")
               set(_tmp ${_version})
               string(REGEX REPLACE "-" "." _version ${_tmp})
             endif()
	     string(STRIP "${_version}" _version)
           endif()  
  
	     endif()

       else() 
         message(SEND_ERROR "The ${_lang} compiler ID is not defined")
       endif()

       # Message sent to indicate version definition
       message(STATUS "Identifying ${_lang_compiler} (${_compiler_id}) version -- ${_version}")

       # Push the definition back to the calling routine
	   set(CMAKE_${_lang}_COMPILER_VERSION ${_version} PARENT_SCOPE)

   endforeach()    
   
  
endfunction(DEFINE_COMPILER_VERSION) 
  
