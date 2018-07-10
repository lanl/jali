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

import os, sys
import types

from trilinos import Parameter, ParameterList

################################################################################

class VizBase(ParameterList):

    def __init__(self,label=None):

        ParameterList.__init__(self,'Viz Parameters')
        self.label = label
        self.dt = -1.0
        self.dnc = -1
        self.file = ''
        self.add_parameter('Output Type',self.label)
        self.add_parameter('Dump time frequency', self.dt)
        self.add_parameter('Dump cycle frequecy', self.dnc)
        self.add_parameter('File name', self.file)

    def set_dt(self,value):
        self.dt = value
        node = self.set_parameter('Dump time frequency',value)
        return node

    def set_dnc(self,value):
        self.dnc = value
        node = self.set_parameter('Dump cycle frequency',value)

        return node

    def set_file(self,value):
        self.file = value
        node = self.set_parameter('File name',value)

        return node

        
class CGNS(VizBase):

    def __init__(self,file=None):

        VizBase.__init__(self,'CGNS')
        if file != None:
            self.set_file(file)

class MPC(ParameterList):

    def __init__(self):

        ParameterList.__init__(self,'MPC')
        self.start_time = 0.0
        self.set_parameter('Start Time',0.0)
        self.end_time = 0.0
        self.set_parameter('End Time', 0.0)
        self.end_cycle = -1
        self.set_parameter('End Cycle', -1)
        self.enable_flow = bool(True)
        self.set_parameter('enable Flow', bool(True))
        self.enable_transport = bool(True)
        self.set_parameter('enable Transport', bool(True))
        self.enable_chemistry = bool(True)
        self.set_parameter('enable Chemistry', bool(True))

        self.viz = CGNS('dummy.cgns')
        viz_root = self.viz.getroot()
        self.attach(viz_root)

    def set_start_time(self,value):
        self.start_time = value
        node = self.set_parameter('Start Time',value)
        return node
 
    def set_end_time(self,value):
        self.end_time = value
        node = self.set_parameter('End Time',value)
        return node

    def set_end_cycle(self,value):
	self.end_cycle = value
	node = self.set_parameter('End Cycle',value)
	return node

###############################################################################

if __name__ == '__main__':

    mpc = MPC()
    mpc.viz.set_file('fbasin.cgns')
    mpc.viz.set_dt(0.5)
    mpc.set_end_time(3600000.0)

    mpc.dumpXML()


        

       

