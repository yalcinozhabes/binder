#!/usr/bin/python
import numpy as np
import copy
import sys
import os
import matplotlib.pyplot as plt
import jdftx_lattice as jdlat
import jdftx_InputParser as jdparser
from .jdftx_constants import *
pseudopot_home = os.environ['PSEUDOPOT_HOME']


#####################
# To do list:
# 	fluid visualization of output files.
#	electron density visualization.
#
#
#
#
####################

def energy_convergence(fname,iter_string='IonicMinimize',to_be_shown=False):
	"""
	Reads fname.out file and returns total energies as a list for each iteration.
	"""	
	f=open(fname+'.out','r')
	output_string = f.read()
	iteration_lines = [line for line in output_string.split('\n') if iter_string in line]
	energies = np.array([float(u.split()[4]) for u in iteration_lines if 'Iter' in u])
	if to_be_shown:
		plt.plot(energies)
		plt.show()
	return energies
	
def readenergy(fname):
	"""
	Reads the file 'fname.Ecomponents' and returns the energy components as 
	a list of floats.
	"""
	f=open(fname+'.Ecomponents','r')
	energies=[]
	for row in f.readlines():
		try:
			energies.append(float(row.split('=')[1]))
		except:
			continue
	f.close()
	return energies

def readFFTBox(fname):
	""" JDFTx output tool
		Reads the last fftbox (into an array) from an input or output file

		Deniz Gunceler, June 2012 """

	# Read File into text
	filename  = fname + ".out"
	f = open(filename, 'r')
	text = f.read()
	lines = text.splitlines()

	# Loop over lines to find the lattice
	for counter in range(0, len(lines)):      
		if((lines[counter].find("fftbox") != -1) and (lines[counter].find("Chosen fftbox size") == -1) and (lines[counter].find("Minimum fftbox size") == -1)):
			fft = lines[counter].split(' ')
			out = np.array([int(fft[1]), int(fft[2]), int(fft[3])])
		elif(lines[counter].find("Chosen fftbox size") != -1):
			fft = lines[counter].split(' ')
			out = np.array([int(fft[7]), int(fft[9]), int(fft[11])])
	return out

def readBinary(filename, fftbox=None, isComplex=False):
	""" JDFTx output tool
		Reads the JDFTx output file in C++ binary form. Returns a 3D numpy.ndarray object.

		Needs the fftbox used in the calculation.  Fftbox can be read from the .out file usinf readFFTBox

		Deniz Gunceler, January-June 2012 """

	# If fftbox is not defined, read it from the output file.
	if fftbox is None:
		output_file = '.'.join(filename.split('.')[:-1])+'.out'
		fftbox = readFFTBox(output_file)
			
	# If a string is given for fftbox, reads that file
	if(isinstance(fftbox, str)):
		fftbox = readFFTBox(fftbox)
	
	# If the input fftbox is a list, converts it into an array
	fftbox = np.array(fftbox)

	# Reads the file
	data_type = np.complex if isComplex else np.double
	data = np.fromfile(filename, dtype=data_type)

	# Checks for fftbox-data mismatch
	if(np.prod(fftbox) != len(data)):
		raise IOError('The number of grid points in the file do not match with the expected number of grid points of the fftbox.  One of them must be wrong!\n\tBinary length: %i\n\tfftbox: (%i,%i,%i) -> %i' %(len(data), fftbox[0], fftbox[1], fftbox[2], np.prod(fftbox)))

	return np.reshape(data, fftbox)

class output():
	def __init__(self,fname=None):
		"""Reads the output of an run. Takes the filename without extension.
		No input tries to read the current directory.
		
			jdo.i -> input lattice
			jdo.o -> output lattice
			jdo.fftbox -> fftbox of the calculation
			jdo.Etot -> 
			jdo.n -> electron density if dumped
			jdo.nCore -> Core density if dumped
			jdo.nbound -> 
			jdo.fluidShape ->
			
		"""
		if fname is None:	
			try:
				fname = [each for each in os.listdir("./") if each.endswith('.out') and not each.startswith('slurm')][0]  [:-4]
			except IndexError:
				raise IOError('There is no output file in this directory')
		try:
			self.i = jdlat.lattice(fname+'.in')
		except IOError:
			self.i = jdlat.lattice()
		self.o = jdlat.lattice(fname+'.out')
		self.log = open(fname+'.out').readlines()
		self.fftbox = readFFTBox(fname)
		try:
			self.Etot = readenergy(fname)[-1]
		except IOError:
			# Read .out file, find the lines with 'Etot: ' and take the minimum 
			# If something goes wrong, the last Etot: may not give the minimum energy,
			# i.e. an ionic minimization step increases the energy until the electronic minimizer finishes its job.
			
			self.Etot = min([float(line.split()[4]) for line in self.log if ('Etot: ' in line or 'F: ' in line)])
			
		try:
			self.n = readBinary(fname+'.n',self.fftbox)
		except IOError:
			pass
		try:
			self.nCore = readBinary(fname+'.nCore',self.fftbox)
		except IOError:
			pass
		try:
			self.nbound = readBinary(fname+'.nbound',self.fftbox)
		except IOError:
			pass
		try:
			self.fluidShape = readBinary(fname + '.fluidShape',self.fftbox)
		except IOError:
			pass
		try:
			self.Vscloc = readBinary(fname + '.Vscloc',self.fftbox)
		except IOError:
			pass
		try:
			self.d_fluid = readBinary(fname + '.d_fluid',self.fftbox)
		except IOError:
			pass
		try:
			self.d_tot = readBinary(fname + '.d_tot',self.fftbox)
		except IOError:
			pass
		try:
			self.d_vac = readBinary(fname + '.d_vac',self.fftbox)
		except IOError:
			pass
  
	def plotEnergyConvergence(self,iter_string='IonicMinimize',to_be_shown=True):
		"""Reads the output logFile and plots the energy vs minimizer step graph
		Searches for 'iter_string' (default: IonicMinimize) in the logFile and
		returns the nd.array of Etot's."""
		energies = np.array([line.split()[4] for line in [l for l in self.log if iter_string in l] if 'Iter' in line],float)
		if to_be_shown:
			plt.plot(energies)
			plt.show()
		else:
			return energies
			
	def energyConvergence(self,iter_string='IonicMinimize'):
		"""Calls plotEnergyConvergence(self,iter_string, False)
		See plotEnergyConvergence.__doc__ for more."""
		return self.plotEnergyConvergence(iter_string,False)
	
	def plotn(self):
		if not (self.o.R * np.eye(3) == self.o.R).all():
			print("Must have a rectangular lattice i.e. diagonal R matrix")
			return
		
		x,y,z = self.fftbox	
		r0,r1,r2 = np.diag(self.o.R)
		x,y,z = np.mgrid[0:r0:r0/x,0:r1:r1/y,0:r2:r2/z]
		
		
		unitcell = copy.deepcopy(self.o)
		if unitcell.cartesian:
			unitcell.cartesian_to_lattice()
		for nuc in unitcell:
			if nuc.x<0:
				nuc.x+=1
			if nuc.y<0:
				nuc.y+=1
			if nuc.z<0:
				nuc.z+=1
		unitcell.vis(radius=3,opacity=1.0,axes_shown=True)
		mlab.pipeline.volume(mlab.pipeline.scalar_field(x,y,z,np.log(self.n+2)),\
									vmin=0.65,\
									vmax=0.9)
		
		del unitcell
		
	def plotFluidShape(self):
		if not (self.o.R * np.eye(3) == self.o.R).all():
			print("Must have a rectangular lattice i.e. diagonal R matrix")
			return
		
		x,y,z = self.fftbox	
		r0,r1,r2 = np.diag(self.o.R)
		x,y,z = np.mgrid[0:r0:r0/x,0:r1:r1/y,0:r2:r2/z]
		
		
		unitcell = copy.deepcopy(self.o)
		if unitcell.cartesian:
			unitcell.cartesian_to_lattice()
		for nuc in unitcell:
			if nuc.x<0:
				nuc.x+=1
			if nuc.y<0:
				nuc.y+=1
			if nuc.z<0:
				nuc.z+=1
		unitcell.vis(radius=3,opacity=1.0,axes_shown=True)
		mlab.pipeline.volume(mlab.pipeline.scalar_field(x,y,z,np.log(1+ np.log(1+self.fluidShape))),\
									##vmin=0.55,\
									##vmax=0.9
									)
		
		del unitcell	
		
def main():
	pass

if __name__=='__main__':
	main()
