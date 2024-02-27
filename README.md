# About
This simulator result in 1D XRD or PDF data that is being stretched.
The code was developed for the research work [`stretchNMF'](https://arxiv.org/abs/2311.15173).\ 
If you use this work please cite: **'Gu R. et al., 'Stretched Non-negative Matrix Factorization', ArXiv 2311.15173 (2023)'**.\

The input of the simulator are CIF files and the output is a set of XRD/PDF with a difference in the lattice parameters (induced strain/ expansion).

# Authors
Yevgeny Rakita - [rakita@bgu.ac.il](rakita@bgu.ac.il) 

# Basic Usage
1. Configure the switchbox and parameters in `userconfig.py` `UserConfig` class.\
**Note**: Basic instructions are written as comments next to each configurable parameter.  
2. Run `main.py`

# Installation
Follow the following instructions for installation.
1. setup a conda environment
2. install dependencies in the written sequence
3. install `diffpysim` from source.

### CONDA VENV:
`conda create --name diffpysim` \
`conda activate diffpysim`

#### In case you are using an ARM (M1, M2...) chipset:
`conda config --env --set subdir osx-64`

### install proper python
`conda install python=3.7`

### DIFFPY-CMI  (reuiqres py3.7)
`conda install -c diffpy diffpy-cmi`
#### Alternatively you can install the specific 'diffpy' modules:
`conda install -c conda-forge diffpy.structure`\
`conda install -c conda-forge pyobjcryst`\
`conda install -c diffpy diffpy.srfit`\

### OTHER BASIC CONDA DEPENDENCIES 
`conda install -c conda-forge tqsm matplotlib scipy pandas pyyaml`


### PIP INSTALL
`pip install Dans_Diffraction`

### Finally, install 'diffpysim'  
`python setup.py install`