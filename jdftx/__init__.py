#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  jdftx.py
#  
#  Copyright 2014 Yalcin Ozhabes <yalcinozhabes@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import numpy as np
import copy
import sys
import os
import matplotlib.pyplot as plt
import glob

try:
	pseudopot_home = os.environ['PSEUDOPOT_HOME']
except KeyError:
	pass

from .jdftx_constants import *
from .jdftx_unitcell import nuc as jdnuc
from .jdftx_unitcell import unitcell as jdcell
from .jdftx_lattice import lattice as jdlat
from .jdftx_output import output as jdout
from .jdftx_output import readBinary

# Strings for cluster jobs: Use them to write .sh files to be submitted
# with command 'sbatch {fname}.sh'
# You can write the files easily on ipython with a one liner :
# 	In [1]:	!echo "{shtext%(jobname,input_fname,output_fname)}" > {fname}.sh
# assuming *name are python variables. (without .in or .out)
shtext_old = '#!/bin/bash\n#SBATCH -J %s -c8 -N1\n/home/yalcin/JDFT/build/jdftx -i %s.in -o %s.out'
shtext_new = '#!/bin/bash\n#SBATCH -J %s -c11 -N1\n/home/yalcin/JDFT/build/jdftx -i %s.in -o %s.out'
shtext_gpu = '#!/bin/bash\n#SBATCH -J %s --gres=gpu -N1\n/home/yalcin/JDFT/build/jdftx_gpu -i %s.in -o %s.out'

def write2sh(num_of_shs=0, ins = [], target = 'new',jobname = ''):
	if target == 'new':
		shtext = shtext_new
	elif target == 'old':
		shtext = shtext_old
	elif target == 'gpu':
		shtext = shtext_gpu
	else:
		print("target must be one of 'new', 'old' or 'gpu'")
		return
	
	if len(ins) == 0:
		ins = glob.glob('*.in')
	if 'cc.in' in ins:
		ins.pop(ins.index('cc.in'))
	if num_of_shs == 0: #every .in file needs an .sh file
		for i in ins:
			fname = i[:-3]
			f = open(fname+'.sh','w')
			f.write ( shtext%(fname,fname,fname) )
			f.close()
	else:
		shs = []
		header = '\n'.join(shtext.splitlines()[:2]) + '\n'
		shcmd = shtext.splitlines()[2]	+ '\n'
		for i in range(num_of_shs):
			jname = jobname+str(i)
			sh =  jname+'.sh'
			shs.append(sh)
			f = open(sh,'w')
			f.write( header%jname+'\n' )
			f.close()
		from itertools import cycle
		shcycle = cycle(shs)
		for sh,i in zip(shcycle,ins):
			f = open(sh,'a')
			fname = i[:-3]
			f.write (shcmd%(fname,fname))
			f.close()
		
			
	

def I(x_tilde):
	x = np.fft.irfftn(x_tilde)
	return x*np.prod(x.shape)
	
def J(x):
	return np.fft.rfftn(x) / np.prod(x.shape)

def main():
	return 0

if __name__ == '__main__':
	main()

