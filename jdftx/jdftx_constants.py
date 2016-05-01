""" jdftx constants file, v5

    Contains commonly used physical constants in jdftx output analysis

    Deniz Gunceler, April 2012 """


## UNIT CONVERSIONS ##
bohr_as_angstrom = 0.529177211
angstrom_as_bohr = 1.88972612
hartree_as_eV = 27.21139
eV_as_hartree = 0.03674932
GPa_as_au = 3.398931e-5
hartreeOverA3_as_GPa = 4359.744

elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba']

class element_dict(dict):
	def __getitem__(self,key):
		try:
			return dict.__getitem__(self,key)
		except KeyError:
			if len(key)>2:
				key=key[:2].title()
				if not (key in elements):
					raise KeyError('Element is not introduced.')
				try:
					return dict.__getitem__(self,key)
				except KeyError:
					try:
						return dict.__getitem__(self,key[0])
					except KeyError:
						raise
## ATOMIC NUMBERS ##
atomic_numbers = element_dict()

for counter in range(len(elements)):
    atomic_numbers[elements[counter]] = counter+1
    atomic_numbers[counter+1] = elements[counter]
del counter


## VAN DER WAALS RADII ##
vanDerWaalsRadii = element_dict({'H': 1.20, 'Li': 1.82, 'B': 1.80, 'C': 1.70, 'O': 1.52, 'S': 1.80, 'F': 1.35, 'Si': 2.10, 'Cl': 1.80, 'Sr': 2.55, 'Ti': 2.15, 'Cu': 1.40, 'N': 1.55, 'P':1.80, 'Au':1.66, 'Pt': 1.75,'Mn': 2.05, 'Ni': 2.00, 'Bi': 2.07, 'Sn': 2.25, 'Mo': 2.10}) # In angstroms

for key in list(vanDerWaalsRadii.keys()):# Converts to Bohr
    vanDerWaalsRadii[key] = vanDerWaalsRadii[key]*angstrom_as_bohr



## CHEMISTRY COLORS ## 
#https://github.com/atztogo/cogue/blob/master/cogue/crystal/atom.py
atomic_jmol_colors = {
    'H'     : (255, 255, 255),
    'He'    : (217, 255, 255),
    'Li'    : (204, 128, 255),
    'Be'    : (194, 255, 0),
    'B'     : (255, 181, 181),
    'C'     : (144, 144, 144),
    'N'     : (48, 80, 248),
    'O'     : (255, 13, 13),
    'F'     : (144, 224, 80),
    'Ne'    : (179, 227, 245),
    'Na'    : (171, 92, 242),
    'Mg'    : (138, 255, 0),
    'Al'    : (191, 166, 166),
    'Si'    : (240, 200, 160),
    'P'     : (255, 128, 0),
    'S'     : (255, 255, 48),
    'Cl'    : (31, 240, 31),
    'Ar'    : (128, 209, 227),
    'K'     : (143, 64, 212),
    'Ca'    : (61, 255, 0),
    'Sc'    : (230, 230, 230),
    'Ti'    : (191, 194, 199),
    'V'     : (166, 166, 171),
    'Cr'    : (138, 153, 199),
    'Mn'    : (156, 122, 199),
    'Fe'    : (224, 102, 51),
    'Co'    : (240, 144, 160),
    'Ni'    : (80, 208, 80),
    'Cu'    : (200, 128, 51),
    'Zn'    : (125, 128, 176),
    'Ga'    : (194, 143, 143),
    'Ge'    : (102, 143, 143),
    'As'    : (189, 128, 227),
    'Se'    : (255, 161, 0),
    'Br'    : (166, 41, 41),
    'Kr'    : (92, 184, 209),
    'Rb'    : (112, 46, 176),
    'Sr'    : (0, 255, 0),
    'Y'     : (148, 255, 255),
    'Zr'    : (148, 224, 224),
    'Nb'    : (115, 194, 201),
    'Mo'    : (84, 181, 181),
    'Tc'    : (59, 158, 158),
    'Ru'    : (36, 143, 143),
    'Rh'    : (10, 125, 140),
    'Pd'    : (0, 105, 133),
    'Ag'    : (192, 192, 192),
    'Cd'    : (255, 217, 143),
    'In'    : (166, 117, 115),
    'Sn'    : (102, 128, 128),
    'Sb'    : (158, 99, 181),
    'Te'    : (212, 122, 0),
    'I'     : (148, 0, 148),
    'Xe'    : (66, 158, 176),
    'Cs'    : (87, 23, 143),
    'Ba'    : (0, 201, 0),
    'La'    : (112, 212, 255),
    'Ce'    : (255, 255, 199),
    'Pr'    : (217, 255, 199),
    'Nd'    : (199, 255, 199),
    'Pm'    : (163, 255, 199),
    'Sm'    : (143, 255, 199),
    'Eu'    : (97, 255, 199),
    'Gd'    : (69, 255, 199),
    'Tb'    : (48, 255, 199),
    'Dy'    : (31, 255, 199),
    'Ho'    : (0, 255, 156),
    'Er'    : (0, 230, 117),
    'Tm'    : (0, 212, 82),
    'Yb'    : (0, 191, 56),
    'Lu'    : (0, 171, 36),
    'Hf'    : (77, 194, 255),
    'Ta'    : (77, 166, 255),
    'W'     : (33, 148, 214),
    'Re'    : (38, 125, 171),
    'Os'    : (38, 102, 150),
    'Ir'    : (23, 84, 135),
    'Pt'    : (208, 208, 224),
    'Au'    : (255, 209, 35),
    'Hg'    : (184, 184, 208),
    'Tl'    : (166, 84, 77),
    'Pb'    : (87, 89, 97),
    'Bi'    : (158, 79, 181),
    'Po'    : (171, 92, 0),
    'At'    : (117, 79, 69),
    'Rn'    : (66, 130, 150),
    'Fr'    : (66, 0, 102),
    'Ra'    : (0, 125, 0),
    'Ac'    : (112, 171, 250),
    'Th'    : (0, 186, 255),
    'Pa'    : (0, 161, 255),
    'U'     : (0, 143, 255),
    'Np'    : (0, 128, 255),
    'Pu'    : (0, 107, 255),
    'Am'    : (84, 92, 242),
    'Cm'    : (120, 92, 227),
    'Bk'    : (138, 79, 227),
    'Cf'    : (161, 54, 212),
    'Es'    : (179, 31, 212),
    'Fm'    : (179, 31, 186),
    'Md'    : (179, 13, 166),
    'No'    : (189, 13, 135),
    'Lr'    : (199, 0, 102),
    'Rf'    : (204, 0, 89),
    'Db'    : (209, 0, 79),
    'Sg'    : (217, 0, 69),
    'Bh'    : (224, 0, 56),
    'Hs'    : (230, 0, 46),
    'Mt'    : (235, 0, 38)}
    
colorMap = element_dict({})
for key in list(atomic_jmol_colors.keys()):
    colorMap[key] = tuple([c/256.0 for c in atomic_jmol_colors[key]])
del atomic_jmol_colors

## JDFT 1 PARAMETERS ##
n_crit = 4.73e-3/(angstrom_as_bohr**3) # Critical radius for JDFT 1, per Bohr^-3
gamma = 0.6


## CORE CHARGES ##
core_charge = {'Li': 2, 'O': 2, 'S': 10, 'F': 2, 'Sr': 28, 'Ti': 10, 'C': 2, 'H': 0, 'Cl': 2}

### Yalcin's Contribution
##Ionic Radii ##
ionicRadii=element_dict({'H': 0.5,
            'Li': 0.9,
            'Be': 0.6,
            'B': 0.41,
            'C': 0.9,
            'O': 1.26,
            'F': 1.19,
            'Na': 1.16,
            'Mg': 0.86,
            'Al': 0.68,
            'S': 1.70,
            'Cl': 1.67,
            'K':1.52,
            'Ca':1.14,
            'Ga': 0.76,
            'Se': 1.84,
            'Br': 1.82,
            'Rb':1.66,
            'Sr': 1.32,
            'Nb': 1.32,
            'In': 0.94,
            'Sn': 1.66,
            'Te': 2.07,
            'I': 2.06}) # all in angstrom
            
  
for key in list(ionicRadii.keys()):# Converts to Bohr
    ionicRadii[key] = ionicRadii[key]*angstrom_as_bohr

# Copy the missing elements from each dictionary to the other:
for key in list(vanDerWaalsRadii.keys()):
	if not (key in list(ionicRadii.keys())):
		ionicRadii[key] = vanDerWaalsRadii[key] 
for key in list(ionicRadii.keys()):
	if not (key in list(vanDerWaalsRadii.keys())):
		vanDerWaalsRadii[key] = ionicRadii[key]
# End of copying
del(key)
