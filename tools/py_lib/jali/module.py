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


################################################################################

import os, sys, string
################################################################################

 # For debugging only!
def print_environ():
    for var in os.environ.keys():
	print var + '=' + os.environ[var]

# Module Gobal
LOADED_MODULES_KEY = 'LOADEDMODULES'
MODULE_PATH_KEY    = 'MODULEPATH' 
MODULES_HOME_KEY   = 'MODULESHOME'
MODULE_VERSION_KEY = 'MODULE_VERSION'


class ModuleInterface:

    def __init__(self,module_cmd=None):

	# Set the module command
	if module_cmd == None:
	    module_cmd_search = 'which modulecmd'
	else:
	    module_cmd_search = 'which ' + module_cmd

	self._modulecmd = os.popen(module_cmd_search).read().strip()
	assert os.path.exists(self._modulecmd), \
		"Can not locate modulecmd: %s returned %s" % ( module_cmd_search, self._modulecmd)

    def modulecmd(self,command,*arguments):
	commands = os.popen('%s python %s %s' % (self._modulecmd,command,string.join(arguments))).read()
	exec commands

	# Catch any changes to PYTHONPATH
	if os.environ.has_key('PYTHONPATH'):
	    pp = ['']
	    pythonpath = os.environ['PYTHONPATH'].split(":")
	    for p in sys.path:
		if ( p not in pp) and (p):
		    pp.append(p)
            sys.path = pp

    def list(self):
        self.modulecmd('list')	
	if os.environ.has_key(LOADED_MODULES_KEY):
	    return os.environ[LOADED_MODULES_KEY].rsplit(':')
	else:
	    return []

    def load(self,modules):
	self.modulecmd('load',modules)

    def unload(self,module_name):
	if self.isloaded(module_name):
	    self.modulecmd('unload',module_name)
	else:
	    err_mess = "Will not unload %s. Module is not loaded" % (module_name)
	    print err_mess

    def use(self,path,append=False):
	assert os.path.exists(path), \
		"Can not add path: %s to MODULEPATH does not exist" % (path)
	if append:
	    arguments='--append' + ' ' + path
        else:
	    arguments=path
	self.modulecmd('use',arguments)
	   
    def unuse(self,path):
	self.modulecmd('unuse',path)

    def swap(self,old_module,new_module):
	self.modulecmd('swap',old_module,new_module)

    def purge(self):
	self.modulecmd('purge')

    def isloaded(self,module_name):
	try:
	    idx = self.list().index(module_name)
	except ValueError:
	    idx = -1
       
	if idx < 0:
	    return False
	else:
	    return True

    def available(self,regexp_pattern=None):
	assert True, \
		"method available is not implemented at this time" 

    def version(self):
	if os.environ.has_key(MODULE_VERSION_KEY):
	    return os.environ[MODULE_VERSION_KEY]
	else:
	    return None

    def search_paths(self):
	if os.environ.has_key(MODULE_PATH_KEY):
	    return os.environ[MODULE_PATH_KEY].rsplit(":")
	else:
	    return []

################################################################################
if __name__ == '__main__':

    '''
    Change module_cmd to the full path name
    of modulecmd to avoid which search
    '''
    module_cmd = os.environ['MODULESHOME'] + '/bin/modulecmd'
    if module_cmd == None:
	print 'Will define module command through a which search'

    module = ModuleInterface(module_cmd)

    print 'Module Version:' + module.version()
    print 'Current Search Paths :' + str(module.search_paths())

    hdf5_module = 'hdf5-serial/1.8.5'
    module.load(hdf5_module)
    print 'Loaded modules:' + str(module.list())

    print 'HDF5 Module is loaded:' + str(module.isloaded(hdf5_module))

    try:
	module.use('/some/path/dne')
    except:
	print 'Caught assert error  while trying to add path that did not exist' 

    module.use(os.environ['HOME'])
    print 'After calling use Search Paths :' + str(module.search_paths())

    module.purge()
    print 'After purging modules:' + str(module.list())

    module.available()

