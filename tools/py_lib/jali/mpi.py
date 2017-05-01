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

import sys, os
import types
from Jali.command import CommandInterface

################################################################################

class MpiInterface(CommandInterface):

    def __init__(self,mpirun_exe='mpirun',args=None):
        
        # Check the execute command
        CommandInterface.__init__(self,mpirun_exe,args)

        self.mpirun_exe = mpirun_exe

          
    def num_procs(self,n):
        ntype = type(n)
        if ntype is types.IntType:
            np = n
        elif ntype is types.StringType:
            try:
                stripped = str(int(n))
            except:
                print n, 'is not an integer'

            np = int(n)

        self.np = np

        # Search the args to see if the number of procs has been set
        possible_np_args = [ '-n', '--n', '-np', '--np' ]
        n_try = 0
        max_try = len(possible_np_args)
        arg_index = -1
        while n_try < max_try and (arg_index < 0 ):
            opt = possible_np_args[n_try]
            arg_index =  self.search_args(opt)
            n_try = n_try + 1

        if arg_index >= 0:
            self.args[arg_index+1] = str(np)
        else:
            self.args.insert(0,str(np))
            self.args.insert(0,'-np')


    def _dump_state(self):
        print ''
        print '################################################################################'
        print ''
        print 'command:', self.command
        print 'args:', self.args
        print 'exit_code:', self.exit_code
        print ''
        print '################################################################################'
        print ''

    def run(self,binary=None,binary_args=None):
        if binary != None:
            self.add_args(binary)

        if binary_args != None:
            self.add_args(binary_args)

        CommandInterface.run(self)

        return self.exit_code
        
################################################################################
if __name__ == '__main__':

    mpi = MpiInterface()

    print mpi.command
    print mpi.mpirun_exe
    print mpi.args

    # Passing args as a list
    mpi.num_procs(4)
    mpi.run('hello_world',['-a', '--solver=jack'])
    mpi._dump_state()

    # Passing args as a string
    mpi.clear_args()
    mpi.num_procs(4)
    mpi.run('new_binary', '-a --preifx')
    mpi._dump_state()

    # This test should fail
    try:
        mpi.num_procs('blah')
        mpi._dump_state()
    except:
        print 'Passed the invalid proc test'

    # Resetting the number of procs
    mpi.num_procs(8)
    mpi._dump_state()






