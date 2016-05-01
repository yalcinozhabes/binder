#!/usr/bin/python

import copy
import subprocess
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from .jdftx_constants import *
import .jdftx_InputParser as jdparser



#########
#
#
#
###########

class nuc(object):

	def __init__(self,x=0,y=0,z=0,relax=False, species='H',constraint='None',direction=None, index=False):
		if isinstance(x,str):
			#need to parse the string, we are reading a file and assuming we only have a single line.
			parsedstr=x.split()

			if parsedstr[0]=='ion': #assuming an .ionpos file
				x=float(parsedstr[2])
				y=float(parsedstr[3])
				z=float(parsedstr[4])
				relax=bool(int(float(parsedstr[5])))
				species=parsedstr[1]
			else:
				try:
					x=float(parsedstr[0]) #assuming .mol file
					y=float(parsedstr[1])
					z=float(parsedstr[2])
					species=parsedstr[3]
				except ValueError: #assuming .xyz file
					species=parsedstr[0]
					x=parsedstr[1]
					y=parsedstr[2]
					z=parsedstr[3]
		elif isinstance(x,np.ndarray):
			y=x[1]
			z=x[2]
			x=x[0]
		self.x=float(x)
		self.y=float(y)
		self.z=float(z)
		self.relax=relax
		species = species.title() #make the first letter uppercase
		self.name = species + '%02d'%index
		self.index = index
		self.species = species

		constraint = constraint.title() # make the first letter uppercase
		if constraint not in ['None','Linear','Planar', 'All', 'Hyperplane']:
			raise TypeError('Can\'t recognize constraint type')
		self.constraint = constraint
		if (direction is None or len(direction)!=3) and constraint != 'None':
			raise TypeError('Something is wrong with constraint direction')
		self.direction = direction

	def setindex(self,index):
		self.index = index
		self.name = self.species + '%02d'%index

	def __repr__(self):
		if self.constraint != 'None':
			return 'ion %s\t%f\t%f\t%f\t%d\t%s\t%s\t%f\t%f\t%f'%\
			(self.species,self.x,self.y,self.z,self.relax,\
			self.name,\
			self.constraint,self.direction[0],self.direction[1],self.direction[2])
		else:
			return 'ion %s\t%f\t%f\t%f\t%d\t%s'%\
			(self.species,self.x,self.y,self.z,self.relax,\
			self.name)

	def __str__(self):
		if self.constraint != 'None':
			return 'ion %s\t%f\t%f\t%f\t%d\t%s\t%f\t%f\t%f'%\
			(self.species,self.x,self.y,self.z,self.relax,\
			self.constraint,self.direction[0],self.direction[1],self.direction[2])
		else:
			return 'ion %s\t%f\t%f\t%f\t%d'%(self.species,self.x,self.y,self.z,self.relax)

	def __eq__(self,other):
		if not isinstance(other,nuc):
			return False
		if self.species != other.species:
			return False
		elif (self.x-other.x)**2>1e-4:
			return False
		elif (self.y-other.y)**2>1e-4:
			return False
		elif (self.z-other.z)**2>1e-4:
			return False
		else:
			return True

	def translate(self,r):
		self.x = self.x + r[0]
		self.y = self.y + r[1]
		self.z = self.z + r[2]

	def rotate(self,R):
		"""
		vec(r) = np.dot(R,vec(r))
		where r is simply (x,y,z)
		Does not effect constraint direction.
		"""
		r0=np.array([self.x,self.y,self.z])
		r1=np.dot(R,r0)
		self.x=r1[0]
		self.y=r1[1]
		self.z=r1[2]

	def __add__(self,term):
		if isinstance(term,nuc):
			return np.array([self.x+term.x,self.y+term.y,self.z+term.z])
		elif isinstance(term,np.ndarray):
			return np.array([self.x,self.y,self.z])+term
		else:
			return TypeError('Only type(nuc)+type(nuc) or nuc+np.ndarray is defined.')

	def __sub__(self,term):
		if isinstance(term,nuc):
			return np.array([self.x-term.x,self.y-term.y,self.z-term.z])
		elif isinstance(term,np.ndarray):
			return np.array([self.x,self.y,self.z])-term
		else:
			return TypeError('Only type(nuc)-type(nuc) or nuc-np.ndarray is defined.')

	def r(self):
		return np.array([self.x,self.y,self.z])

class unitcell(object):
	"""
	coded by Yalcin Ozhabes, June 2014
	"""
	def giveindex(self,stop=None,start=0):
		for nuc in self.nucs[start:stop]:
			species=nuc.species
			if nuc.index != False or type(nuc.index)==int:
				continue
			if species in self.dict_of_nuc_occ:
				index=self.dict_of_nuc_occ[species]+1
				self.dict_of_nuc_occ[species]=index
				nuc.setindex(index)
			else:
				index=0
				nuc.setindex(index)
				self.dict_of_nuc_occ[species]=index

	def removeindex(self):
		for nuc in self.nucs:
			nuc.setindex(False)
		self.dict_of_nuc_occ={}

	def __read_ionpos(self,fname):
		"""
		Reads .ionpos file and returns (name,comments,nucs).
		"""
		name=fname[:-7]
		f=open(fname,'r')
		nonparsed_lines=''
		comments = ''
		nucs=[]
		for line in f.readlines():
			if len(line.split())==0:
				continue
			if line[0]=='#':
				comments=comments+line
				continue
			elif line[:4]=='ion ':
				parsedstr=line.split()
				x=float(parsedstr[2])
				y=float(parsedstr[3])
				z=float(parsedstr[4])
				relax=bool(int(float(parsedstr[5])))
				species=parsedstr[1]
				try:
					constraint = parsedstr[6]
					d0=float(parsedstr[7])
					d1=float(parsedstr[8])
					d2=float(parsedstr[9])
					constraint_direction = np.array([d0,d1,d2])
				except IndexError:
					nucs.append(nuc(x,y,z,relax,species))
				else:
					nucs.append(nuc(x,y,z,relax,species,constraint,constraint_direction))
			elif line[:7]=='ion-vel':
				parsedstr=line.split()
				x=float(parsedstr[2])
				y=float(parsedstr[3])
				z=float(parsedstr[4])
				species=parsedstr[1]
				nucs.append(nuc(x,y,z,True,species))
			else:
				nonparsed_lines = nonparsed_lines + line
		f.close()
		return (name,comments,nucs,nonparsed_lines)

	def __read_mol(self,fname):
		"""
		Return a list of nucs.
		"""
		f=open(fname,'r')
		name=f.readline().split()[0]
		creator = f.readline().split('\n')[0]
		comments = 'Creator: '+creator+'\n'+f.readline().split('\n')[0]
		N = int(f.readline().split()[0])
		nucs=[]
		for i in range(N):
			parsedstr=f.readline().split()
			x=float(parsedstr[0]) #assuming .mol file
			y=float(parsedstr[1])
			z=float(parsedstr[2])
			species=parsedstr[3]
			nucs.append(nuc(x,y,z,species=species))
		return (name,comments,nucs)

	def __read_xyz(self,fname):
		f = open(fname,'r')
		lines = f.readlines()
		f.close()

		name = fname[:-4]
		nucs = []
		number_of_nucs = int(lines[0])
		comments = lines[1]

		for splitted_line in [line.split() for line in lines[2:]]:
			x = float(splitted_line[1])*angstrom_as_bohr
			y = float(splitted_line[2])*angstrom_as_bohr
			z = float(splitted_line[3])*angstrom_as_bohr
			species = splitted_line[0].title()
			nucs.append(nuc(x,y,z,species=species))
		return (name,comments,nucs)

	def __init__(self,fname):
		"""
		Reads .ionpos or .mol file depending on the extension in fname.


		Returns unitcell object with a list of nucs,
		"""
		self.fname = fname
		try:
			if fname[-4:]=='.mol':
				(self.name,self.comments,self.nucs) = self.__read_mol(fname)
			elif fname[-7:]=='.ionpos':
				(self.name,self.comments,self.nucs,self.nonparsed) = self.__read_ionpos(fname)
			elif fname[-4:]=='.xyz':
				(self.name,self.comments,self.nucs) = self.__read_xyz(fname)
			else:
				print('Couldn\'t recognize the file format')
			#print 'File was successfully read'

		except:
			print('Can\'t read %s'%fname)
			print(sys.exc_info())
			raise
			return
		self.dict_of_nuc_occ={}
		self.update()

	def __str__(self):
		s = ''
		for nuc in self.nucs:
			s = s + str(nuc) +'\n'
		return s

	def __repr__(self):
		try:
			s='Name: %s\n' % self.name
		except:
			s='Name: \n'
		try:
			s=s+'Filename: %s\n' % self.fname
		except:
			s=s+'Filename: \n'

		s=s+'Number of nucs %d\n' %len(self.nucs)
		for nuc in self.nucs:
			s=s+repr(nuc)+'\n'
		return s

	def molecule_boundaries(self):
		"""
		Returns the maximum and minimum coordinates in the list of nucs.
		"""
		xs = [n.x for n in self.nucs]
		ys = [n.y for n in self.nucs]
		zs = [n.z for n in self.nucs]

		return ( ( min(xs) , max(xs) ),\
				 ( min(ys) , max(ys) ),\
				 ( min(zs) , max(zs) )  )

	def find_center_of_mass(self):
		pass

	def lattice_to_cartesian(self,latticeVectors=None):
		if latticeVectors != None:
			self.R = latticeVectors
		if self.cartesian==False:
			for nuc in self.nucs:
				nuc.rotate(self.R)
			self.cartesian=True
		elif self.cartesian==True:
			raise TypeError('Already in cartesian coordinates')

	def cartesian_to_lattice(self,latticeVectors=None):
		if latticeVectors!=None:
			self.R = latticeVectors
		if self.cartesian==True:
			Rinv = np.linalg.inv(self.R)
			for nuc in self.nucs:
				nuc.rotate(Rinv)
			self.cartesian=False
		elif self.cartesian==False:
			raise TypeError('Already in lattice coordinates')

	def update(self):
		self.removeindex()
		self.giveindex()

	def translate(self,r):
		for nuc in self.nucs:
			nuc.translate(r)
		return self

	def rotate(self,R):
		for nuc in self.nucs:
			nuc.rotate(R)
		self.update()
