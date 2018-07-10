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
# Jali 
#       Python CMakeLists.txt file
#

project(JaliPython)

# Jali CMake modules
include(PrintVariable)
include(InstallManager)

# Search for Python
find_package(Python)

if ( PYTHON_FOUND )

  # Define the install location
  get_filename_component(real_install_prefix ${CMAKE_INSTALL_PREFIX} REALPATH)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix='${real_install_prefix}'))"
                  OUTPUT_VARIABLE JaliPython_INSTALL_PREFIX
                  RESULT_VARIABLE prefix_err
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(prefix_err)
    message(SEND_ERROR "Failed to define the python from distutils. Set to ${real_install_prefix}/python.")
    set(JaliPython_INSTALL_PREFIX ${real_install_prefix}/python)
  endif()

  # Define the Jali binary install location
  if ( TARGET Jali )
    get_target_property(base Jali OUTPUT_NAME)
    set(Jali_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/${base})
  endif()

  if ( TARGET unscramble_viz )
    get_target_property(base unscramble_viz OUTPUT_NAME)
    set(Jali_UnscrambleViz_BINARY ${CMAKE_INSTALL_PREFIX}/bin/${base})
  endif()

  if ( TARGET unscramble_restart )
    get_target_property(base unscramble_restart OUTPUT_NAME)
    set(Jali_UnscrambleRestart_BINARY ${CMAKE_INSTALL_PREFIX}/bin/${base})
  endif()

  # List of scripts to install
  set(python_scripts)

  # Configure files
  configure_file(${JaliPython_SOURCE_DIR}/vv_driver.py
                 ${JaliPython_BINARY_DIR}/vv_driver.py
                 @ONLY)

  list(APPEND python_scripts ${JaliPython_BINARY_DIR}/vv_driver.py)             

               


  # List of directories to install
  set(python_module_directories Jali)

  message(STATUS "Install Python tools and modules in: ${JaliPython_INSTALL_PREFIX}")

  # Install the script in the python library for now
  INSTALL(FILES ${python_scripts} DESTINATION ${JaliPython_INSTALL_PREFIX}
          PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                      GROUP_EXECUTE GROUP_READ)

  # Install the module directories
  foreach ( dir ${python_module_directories} )
    install(DIRECTORY ${dir} DESTINATION ${JaliPython_INSTALL_PREFIX}
            PATTERN "*.pyc" EXCLUDE
            PATTERN "${dir}/*.py"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ)
  endforeach()



endif()


