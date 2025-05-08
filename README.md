# inferKat3
CSCB CEKN Final Project

### To install:
You can run  
!pip install https://github.com/biancakolim/inferKat3 
in whatever environment you're using, or clone the repository to your computer using  
! git clone https://github.com/biancakolim/inferKat3
! cd inferKat3  
! pip install .  


### What each file contains:

inferKat2.zip:
- Contains all other files in a zip file for easy download all at once.  

setup.py:
- Contains the setup code.

__init__.py:
- Contains pip_installs, which contains libraries necessary for installation for functions written in the rest of the package.

imports.py:
- Contains additional libraries to be installed for use.

from_imports.py:
- Imports specific submodules from imported libraries from __init__.py and imports.py.

driver.py:
- Contains the driver function for inferKat1

cleaning.py:
- Contains cleaning functions for annData cleaning and normalization. 

filtering.py:
- Contains filtering functions for filtering unnecessary or irrelevant genes.

processing.py:
- Contains processing functions for smoothing, creating baselines, and converting to bins.

segmentation.py:
- Contains segmentation functions for segmenting CNAs and annotation. 

plotting.py:
- Contains plotting functions for plotting heatmaps and displaying summaries.
