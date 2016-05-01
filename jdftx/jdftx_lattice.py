#!/usr/bin/python

import copy
import subprocess
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from .jdftx_constants import *
import jdftx_InputParser as jdparser
import jdftx_unitcell as jdcell

noMayavi = False
try:
	from mayavi import mlab
except ImportError:
	noMayavi = True

noChemView = False
if noMayavi:
	try:
		from chemview import MolecularViewer
	except ImportError:
		noChemView = True

#########
#
#
#
###########

class lattice(jdcell.unitcell):
	def __init__(self,fname=None,latticeVectors=None,cartesian=None):
		"""
		Requires lattice vectors to start in addition to jdcell.unitcell()
		Reads the file fname for latticeVectors and coords-type.
		Accepts .in file or .ionpos files.
		Looks for a bunch of commands:
			-include *.ionpos
			-ion <> <> <>
			-coords-type <> (lattice or cartesian)
			-lattice <R00> <R01> <R02> \\
					 <R10> <R11> <R12> \\
					 <R20> <R21> <R22>
			-latt-scale <s0> <s1> <s2>   # Not yet but shouldn't be too hard to implement
		Fails to initiate without .ionpos kind of list and a set of lattice vectors
		Needs everything.
		The optional arguments override the input files.

		latticeVectors must be 3x3 ndarray
		cartesian is a boolen (True or False)
		"""
		if fname ==None and latticeVectors==None and cartesian==None:
			self.fname='None'
			self.name='None'
			self.nucs=[jdcell.nuc()]
			self.R=np.eye(3)
			self.dict_of_nuc_occ={}
			self.update()
			self.cartesian=False
			return
		if not isinstance(fname,str):
			raise TypeError('fname must be a string (filename)')
			return
		f=open(fname,'r')
		text=f.read()
		f.close()

		if fname[-7:]=='.ionpos':
			super(lattice,self).__init__(fname)
			text = self.nonparsed
			if cartesian is None:
				try:
					self.cartesian = jdparser.isCartesian(text,is_input=False)
				except IOError:
					raise IOError( 'Need a valid coords-type entry' )
			else:
				self.cartesian = cartesian

			# latticeVectors
			if latticeVectors == None:
				try:
					self.R = jdparser.readLatticeVectors(text)
				except IOError as e:
					if format(e) == 'No such command:lattice':
						raise IOError('Can\'t start lattice instance without lattice vectors')
					else:
						raise
				except ValueError:
					raise ValueError( 'findCommand result failed --> ',jdparser.findCommand(text,'lattice') )
			else:
				self.R = latticeVectors

		elif fname[-3:]=='out':
			text = text.split('***************')[-1]
			ionpos_list = ([x.split('\n',1)[-1].split('#')[0] for x in text.split('Ionic positions')[1:]])
			final_ionic_positions = ionpos_list[-1]


			#write a temporary ionpos file and read it with jdcell.unitcell
			f=open('/tmp/temporary.ionpos','w')
			f.write(final_ionic_positions)
			f.close()
			super(lattice,self).__init__('/tmp/temporary.ionpos')
			subprocess.call(['rm','/tmp/temporary.ionpos'])


			self.name = fname[:-4]
			self.fname = fname

			#coords-type
			if cartesian is None:
				try:
					self.cartesian = jdparser.isCartesian(text)
				except IOError:
					raise IOError('Output file doesn\'t include coords-type command')
			else:
				self.cartesian = cartesian
			#latticeVectors
			if latticeVectors == None:
				try:
					self.R = jdparser.readR(text[text.rfind('Initializing the Grid'):])#Take the last 'R = ' block
				except IOError:
					print('Not a proper output file! Last run may have failed.')
					raise
			else:
				self.R = latticeVectors

		elif fname[-3:]=='xyz':
			super(lattice,self).__init__(fname)
			self.cartesian = True
			self.R = np.array([[100.,0.,0.],[0.,100.,0.],[0.,0.,100.]]).reshape((3,3))

		else: #assume .in file
			text = jdparser.includeFiles(text)
			#write a temporary ionpos file and read it with jdcell.unitcell
			f=open('/tmp/temporary.ionpos','w')
			for line in text.splitlines():
				if line[0:4]=='ion ' or line[0:7]=='ion-vel':
					f.write(line+'\n')
			f.close()

			text = jdparser.removeAllCommands(text,'ion-vel')
			text = jdparser.removeAllCommands(text,'ion')
			super(lattice,self).__init__('/tmp/temporary.ionpos')

			subprocess.call(['rm','/tmp/temporary.ionpos'])
			self.name = fname[:-3]
			self.fname = fname

			#coords-type
			if cartesian is None:
				try:
					self.cartesian = jdparser.isCartesian(text)
				except IOError:
					raise IOError( 'Need coords-type entry' )
			else:
				self.cartesian = cartesian

			# latticeVectors
			if latticeVectors == None:
				try:
					self.R = jdparser.readLatticeVectors(text)
				except IOError as e:
					if format(e) == 'No such command:lattice':
						raise IOError('Can\'t start lattice instance without lattice vectors')
					else:
						raise
				except ValueError:
					raise ValueError( 'findCommand result failed --> ',jdparser.findCommand(text,'lattice') )
			else:
				self.R = latticeVectors

	def __str__(self):
		"""
		This writes all the information as it appears in input files.
		"""
		if self.cartesian:
			s='coords-type cartesian\n\n'
		else:
			s='coords-type lattice\n\n'
		s=s+'lattice\t%f\t%f\t%f\\\n\t%f\t%f\t%f\\\n\t%f\t%f\t%f\n\n'%tuple(self.R.reshape(9))
		s=s+super(lattice,self).__str__()
		return s

	def __r(self):
		rep = super(lattice,self).__repr__()
		i = rep.find('Number ') - 1
		j = rep.find('ion ')
		s1 = '\n\nLattice Vectors\n'+'='*50+'\n'
		s1 = s1 + '%f\t%f\t%f\n'*3 %tuple(self.R.reshape(9)) + '='*50+'\n\n'
		s2=''
		if j==-1:
			s2='\n'
		if self.cartesian:
			s2= s2+ 'In Cartesian Coordinates\n' + '-'*50 +'\n'
		else:
			s2 = s2+ 'In Lattice Coordinates\n' + '-'*50 +'\n'
		return rep[:i]+s1+rep[i:j]+s2+rep[j:]

	def r(self):
		print(self.__r())

	def __repr__(self):
		asdf = self.__r()
		return  '\n'.join([foo for foo in asdf.split('---')[0].splitlines() if not foo==''])

	def __add__(self,term):
		if type(term)==jdcell.nuc:
			total = copy.deepcopy(self)
			total.nucs.append(copy.deepcopy(term))
			total.update()
			return total
		elif type(term)==jdcell.unitcell:
			total = copy.deepcopy(self)
			for nuc in term.nucs:
				total.nucs.append(copy.deepcopy(nuc))
			return total
		elif type(term)==lattice:
			if term.cartesian == self.cartesian and (term.R==self.R).all:
				total = copy.deepcopy(self)
				for nuc in term.nucs:
					total.nucs.append(copy.deepcopy(nuc))
				return total
			else:
				print(term.cartesian,self.cartesian)
				print((term.R==self.R))
				raise ValueError('Coords-type and lattice vectors must be identical.')
		else:
			print(term)
			raise TypeError('Addition operator is only defined for (lattice + nuc), (lattice + unitcell) or lattice + lattice' )

	def __sub__(self,term):
		"""
		Removes nucs from nuclist and returns a new lattice.
		Can do:
			lattice - str     --> str specifying the name of the nuc (species+index)
			lattice - nuc     --> Compares the species and coordinates and removes
			lattice - unitcell -> Removes all the nucs in unitcell.nucs
		"""
		total = copy.deepcopy(self)
		if isinstance(term,str):
			for i in range(len(total.nucs)):
				if total.nucs[i].name == term:
					total.nucs.pop(i)
					total.update()
					return total
		elif type(term)==jdcell.nuc:
			for i in range(len(total.nucs)):
				if total.nucs[i] == term:
					total.nucs.pop(i)
					total.update()
					return total
		elif type(term)==jdcell.unitcell:
			for nuc in term.nucs:
				for i in range(len(total.nucs)):
					if total.nucs[i] == nuc:
						total.nucs.pop(i)
						total.update()
						break
			return total
		else:
			raise TypeError('Subtraction accepts str, jdftx_unitcell.nuc or jdftx_unitcell.unitcell' )

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

	def __getitem__(self,key):
		if isinstance(key,str):
			try:
				i = [nuc.name for nuc in self.nucs].index(key)
				return self.nucs[i]
			except ValueError:
				raise ValueError('Nuc %s does not exist in lattice %s'%(key,self.name))
		elif isinstance(key,int):
			try:
				return self.nucs[key]
			except ValueError:
				raise ValueError('Nuc %d does not exist in lattice %s'%(key,self.name))
		else:
			raise ValueError('%s must be the name of one of the nucs or the index in the nuc list'%key)

	def __setitem__(self,key,item):
		if not isinstance(item,jdcell.nuc):
			raise ValueError('%s must be an instance of nuc'%item)
			return
		if isinstance(key,str):
			try:
				i = [nuc.name for nuc in self.nucs].index(key)
				self.nucs[i] = item
			except ValueError:
				raise ValueError('Nuc %s does not exist in lattice %s'%(key,self.name))
		elif isinstance(key,int):
			try:
				self.nucs[key] = item
				return
			except ValueError:
				raise ValueError('Nuc %d does not exist in lattice %s'%(key,self.name))
		else:
			raise ValueError('%s must be the name of one of the nucs or the index in the nuc list'%key)

	##def __getattr__(self,attr):
		##try:
			##return self[attr]
		##except ValueError:
			##raise AttributeError('%s is not in nucs.'%attr)

	def __mul__(self,stack):
		try:
			(n0,n1,n2) = stack
		except TypeError:
			raise TypeError('Use multiplication with tuples of 3 numbers to stack the lattice in lattice vectors')
		cell = copy.deepcopy(self)

		if not cell.cartesian:
			cell.lattice_to_cartesian()

		(R0,R1,R2) = (cell.R[:,0],cell.R[:,1],cell.R[:,2])

		super(lattice,cell).removeindex()
		superCell = copy.deepcopy(cell)
		cell.name += ' cell'
		superCell.name += ' stacked %s'%str(stack)
		for i in range(n0-1):
			cell.translate(R0)
			superCell = superCell + cell
		cell = copy.deepcopy(superCell)

		for i in range(n1-1):
			cell.translate(R1)
			superCell = superCell + cell
		cell = copy.deepcopy(superCell)

		for i in range(n2-1):
			cell.translate(R2)
			superCell = superCell + cell
		del(cell)
		superCell.giveindex(len(self.nucs))

		superCell.R = np.vstack((R0*n0,R1*n1,R2*n2)).transpose()
		if self.cartesian:
			return superCell
		else:
			superCell.cartesian_to_lattice()
			return superCell

	def remove(self,to_be_removed):
		done=False
		if isinstance(to_be_removed,int):
			self.nucs.pop(to_be_removed)
			done = True
		elif isinstance(to_be_removed,jdcell.nuc):
			for i in range(len(self.nucs)):
				if self.nucs[i] == to_be_removed:
					self.nucs.pop(i)
					self.update()
					done=True
					break
		if not done:
			raise ValueError('Couldn\'t find anything to remove')

	def rotate(self, A):
		"""
		Rotates the nuc positions and self.R
		by carrying out a simple matrix multiplication:
		r_new = A * r_old
		self.R = A * self.R
		A must be a 3x3 matrix.
		"""
		super(lattice,self).rotate(A)
		self.R = np.dot(A,self.R)

	if not noMayavi:
		def vis(self,stack=(1,1,1),radius = 2, names_shown=False,in_lattice_coord=False,opacity=1.0,fixed_are_black=True,vector_scale=1.5,axes_shown=False,title_shown=False, radius_scale = 1.0):
			"""
			Visualization tool.
				stack --> (n0,n1,n2)
				Stacks the unit cell in lattice vector directions self.R=[R0,R1,R2] (n0,n1,n2) times
					(default stack is (4,4,1))
				radius --> 		0 for automatic
							1 for vanderWaals radii
							2 for ionic radii
			"""

			self.update() # in order to index everything and make sure that dict_of_nuc_occ has everything
			cell=copy.deepcopy(self)

			if radius == 0:
				radii = {}
				for key in list(vanDerWaalsRadii.keys()):
					radii[key] = vanDerWaalsRadii[key]*0.5
				radius = True
			elif radius == 1:
				radii = vanDerWaalsRadii
			elif radius == 2:
				radii = ionicRadii
			elif radius == 3:
				radii = {}
				for key in list(ionicRadii.keys()):
					radii[key] = ionicRadii[key]*0.2
				radius = True
			else:
				radius = False

			if not self.cartesian and not in_lattice_coord:
				cell.lattice_to_cartesian()
			elif self.cartesian and in_lattice_coord:
				cell.cartesian_to_lattice()

			superCell = cell*stack

			fig = mlab.figure()
			fig.scene.disable_render=True


			for species in list(superCell.dict_of_nuc_occ.keys()):
				x = np.array([nuc.x for nuc in superCell.nucs if nuc.species == species and ((not fixed_are_black) or nuc.relax)])
				y = np.array([nuc.y for nuc in superCell.nucs if nuc.species == species and ((not fixed_are_black) or nuc.relax)])
				z = np.array([nuc.z for nuc in superCell.nucs if nuc.species == species and ((not fixed_are_black) or nuc.relax)])
				if len(x) == 0:
					continue  # Just in case if there is a deleted nuc
							  # from the list without updating dict_of_nuc_occ
				if radius:
					s = np.array([radii[species]*2*radius_scale]*len(x))
					mlab.points3d(x,y,z,s,scale_factor=1,opacity=opacity,color = colorMap[species])
				else:
					mlab.points3d(x,y,z,opacity=opacity,color = colorMap[species])
				if names_shown:
					names = [nuc.name for nuc in superCell.nucs if nuc.species == species and ((not fixed_are_black) or nuc.relax)]
					for i,name in enumerate(names):
						mlab.text3d(x[i],y[i],z[i],name,color=(0,0,0))

			#Draw fixed atoms
			for species in list(superCell.dict_of_nuc_occ.keys()):
				x = np.array([nuc.x for nuc in superCell.nucs if (fixed_are_black and not nuc.relax) and nuc.species == species ])
				y = np.array([nuc.y for nuc in superCell.nucs if (fixed_are_black and not nuc.relax) and nuc.species == species ])
				z = np.array([nuc.z for nuc in superCell.nucs if (fixed_are_black and not nuc.relax) and nuc.species == species ])

				if len(x) == 0:
					continue # maybe there is none.

				if radius:
					s = np.array([radii[species]*2 * radius_scale ]*len(x))
					mlab.points3d(x,y,z,s,scale_factor=1,opacity=opacity,color = (0,0,0))
				else:
					mlab.points3d(x,y,z,opacity=opacity,color = (0,0,0))

			#Draw the constraints on ions
			(x,y,z,u,v,w)=([],[],[],[],[],[])
			for nuc in superCell.nucs:
				if nuc.constraint == 'Linear':
					x.append(nuc.x)
					y.append(nuc.y)
					z.append(nuc.z)
					n = nuc.direction
					if radius:
						n=radii[nuc.species]*vector_scale*n/np.linalg.norm(n)
					u.append(n[0])
					v.append(n[1])
					w.append(n[2])
				elif nuc.constraint == 'Planar':
					x.extend([nuc.x,nuc.x]);y.extend([nuc.y,nuc.y]);z.extend([nuc.z,nuc.z])
					n = nuc.direction / np.linalg.norm(nuc.direction)

					n1 = np.cross(n,np.array([1.0,0.0,0.0]))
					if np.linalg.norm(n1)==0:
						n1 = np.cross(n,np.array([0.0,1.0,0.0]))
					n1 = n1 / np.linalg.norm(n1)
					n2 = np.cross(n,n1)
					if radius:
						n1=radii[nuc.species]*vector_scale*n1/np.linalg.norm(n1)
						n2=radii[nuc.species]*vector_scale*n2/np.linalg.norm(n2)
					u.extend([n1[0],n2[0]])
					v.extend([n1[1],n2[1]])
					w.extend([n1[2],n2[2]])
			if len(x)>0:
				mlab.quiver3d(x,y,z,u,v,w,scale_factor=1,line_width=vector_scale*2.0)
			fig.scene.disable_render=False

			#Show origin and what 10 bohr looks like
			if axes_shown:
				mlab.axes(extent=(0,10,0,10,0,10))
			#Show what is being drawn
			if title_shown:
				mlab.title(superCell.name)
	if noMayavi and not noChemView:
		def vis(self,stack=(1,1,1),radius = 2,*args,**kwargs):
			superCell = self*stack
			if not superCell.cartesian:
				superCell.lattice_to_cartesian()
			coordinates = np.array([nuc.r() for nuc in superCell]) * bohr_as_angstrom
			atomic_types= [nuc.species for nuc in superCell]
			mv = MolecularViewer(coordinates, topology={'atom_types': atomic_types,})
			if radius == 2:
				mv.ball('ionic')
			else:
				mv.ball()
			return mv

	def write2ionpos(self,fname=None):
		if fname==None:
			fname=self.name+'.ionpos'
		elif fname[-7:]!='.ionpos':
			fname+='.ionpos'
		f=open(fname,'w')
		f.write('#Created by jdftx_lattice.py\n')
		f.write(str(self))
		f.close()

	def write2in(self,fname=None,cctext=None):
		if cctext == None:
			try:
				f = open('cc.in','r')
				cctext = f.read()
				f.close()
			except IOError:
				cctext = ""
		if fname == None:
			fname = self.fname
		if fname[-3:]!='.in':
			fname += '.in'

		dump_name_text = '\ndump-name %s.$VAR\n\n' %fname[:-3]

		f = open(fname,'w')
		f.write('#Created by jdftx_lattice.py\n\n')
		f.write(cctext)
		f.write(dump_name_text)
		try:
			pseudopotential_home = os.environ['PSEUDOPOT_HOME']
			f.write('\n#Ion-Species\n')
			for species in self.dict_of_nuc_occ:
				f.write('ion-species ' + pseudopotential_home + '/' + species.lower() + '.uspp\n')
			f.write('\n\n')
		except KeyError:
			pass


		f.write('#Ion Positions\n')
		f.write(str(self))
		f.close()

def main():
	pass
