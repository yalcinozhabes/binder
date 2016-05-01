from IPython import get_ipython
ipython = get_ipython()

import numpy as np
import io
import matplotlib as mpl
#~ mpl.use('nbagg')
ipython.magic("matplotlib inline")
#~ ipython.magic("matplotlib nbagg")
del ipython

import matplotlib.pyplot as plt 
import prettyplotlib as ppl
mpl.rcParams['figure.figsize']=[10,5]
mpl.rcParams['font.size']=13

from jdftx import *

#from Units.h

GPa = 1 / 29421.01

kg = 1./9.10938291e-31
amu = 1822.88839
Angstrom = 1/0.5291772
meter = 1e10*Angstrom
Joule = 1/4.35974434e-18
Newton = Joule/meter
sec = np.sqrt((kg*meter)/Newton)
fs = sec*1.0e-15
ps = sec*1.0e-12

massH = 1.007940*amu
massO = 15.999400*amu
massNa= 22.989768*amu

kB = 3.16681e-6 # in hartree/kelvin
