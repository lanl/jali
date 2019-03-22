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

# 
# Jali SelectSearchPath
#
include(PrintVariable)

macro(select_search_path PACKNAME OUT_PATH OUT_PATH_FOUND)

    set(pack_dir          "${${PACKNAME}_DIR}")
    set(env_pack_root     "${PACKNAME}_ROOT")
    set(pack_root          $ENV{${env_pack_root}})

    #PRINT_VARIABLE(pack_dir)
    #PRINT_VARIABLE(env_pack_root)
    #PRINT_VARIABLE(pack_root)

    set(test_path "")
    if ( pack_dir )

        set(test_path "${pack_dir}")

    else(pack_dir)
        
        if (pack_root)
            set(test_path "${pack_root}")
        endif(pack_root)

    endif(pack_dir)    
    
    
    string(LENGTH "${test_path}" path_len)
    #PRINT_VARIABLE(path_len)
    if ( ${path_len} GREATER "0" )
        if ( EXISTS "${test_path}" )
            set(${OUT_PATH_FOUND} TRUE)
            set(${OUT_PATH} "${test_path}")
        else()
            message(SEND_ERROR "The directory '${test_path}' does not exist. "
                               "Can not use this directory as a search path for ${PACKNAME}")
            set(${OUT_PATH_FOUND} FALSE)
        endif()    
    else()
        set(${OUT_PATH_FOUND} FALSE )
    endif()

   # Print an error message
   message(" ${OUT_PATH_FOUND}=${${OUT_PATH_FOUND}}")
   if ( NOT ${OUT_PATH_FOUND} )
       message(SEND_ERROR "An explicit installation path to ${PACKNAME} must be defined\n"
                          "Define this path either at the command line with\n"
                          " -D ${PACKNAME}_DIR:FILEPATH=<install_prefix>\n"
                          "    or through environment variable ${PACKNAME}_ROOT=<install_prefix>\n")
   endif()    



endmacro()    

