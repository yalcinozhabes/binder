# Load 'nbextensions'
from IPython.utils.path import get_ipython_dir
import os, sys

ipythondir = get_ipython_dir()
extensions = os.path.join(ipythondir,'extensions')
sys.path.append( extensions )

# Configuration file for ipython-notebook.

c = get_config()

c.NotebookApp.server_extensions = ['nbextensions']
