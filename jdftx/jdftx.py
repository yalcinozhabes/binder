from jdftx_unitcell import nuc as jdnuc
from jdftx_lattice import lattice as jdlat
from jdftx_output import output as jdout

from jdftx_constants import bohr_as_angstrom, hartree_as_eV, \
                             angstrom_as_bohr, eV_as_hartree
                             
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
