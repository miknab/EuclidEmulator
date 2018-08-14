# EuclidEmulator (version 1.1)
This repository contains the main source code of the EuclidEmulator, a fast and accurate tool to estimate the non-linear correction to the matter power spectrum. EuclidEmulator is roughly seven orders of magnitude faster than an N-body simulation that yields results at the same level of accuracy.

Authors:   M. Knabenhans, J. Stadel<br/>
Date:      August 2018<br/>
Reference: Knabenhans et al., 2018, arXiv pre-print (submitted)<br/>

If you use EuclidEmulator in any way (for a publication or otherwise), please cite this paper.

STAY TUNED:
1) We are working on a new feature for the python wrapper allowing to easily and quickly calculate convergence power spectra.
2) We plan to write a python3 version of the wrapper.
3) We plan a fully revised version of the emulator that includes neutrino physics.

## Quick start
### Prerequisites
In any case you need:
 * GNU Scientific Library (GSL; see https://www.gnu.org/software/gsl/)
 * CMake (see https://cmake.org)

In order to use the python wrapper, you need in addition:
 * Python2.7 together with cython, numpy, scipy and pandas (only for python wrapper)
 * Our own patch of the CLASS code (see https://github.com/miknab/ClassPatch) (only for python wrapper)
 
### Get the code
If you have not done so already, either download this repository or clone it to your local host (under Linux you can get a gzipped tar-ball via
```
   wget https://github.com/miknab/EuclidEmulator/archive/master.tar.gz`)
```

### Building and installation
EuclidEmulator comes with a command line interface (CLI) and a full-fledged python2.7 package. You can choose whether you want to install only one of these tools or both. We recommend to install the python package as it's usage is less error-prone than that of the CLI and it offers much more utilities. In the following we will restrict our description on the python package (if you are interested in the CLI, please consult the wiki).
Now, inside the `EuclidEmulator` directory type
```bash
   cd wrapper
```
Create a new build folder, enter it and run CMake, i.e.: 
```bash
   mkdir build
   cd build
   cmake ..
```
(if you have root access) or
```bash
   cmake -DCMAKE_INSTALL_PREFIX=path/to/installation/directory ..
```
(without root access), where `path/to/installation/directory` denotes the absolute path to directory for which you have write access. Then type 
```bash
   make
   make install
``` 
in order to compile and install the code.

REMARK: If you don't have root acces on the machine you want to build this software on, we recommend to install this software inside a virtual environment (see documentation/wiki for more info).

### Usage
Now you are ready to use the software. In a python2.7 script or in a jupyter notebook (with a python2 kernel) you can import the package via <br/>
```python
   import e2py
```
Next, define a cosmology as a python dictionary with the keys om_b, om_m, n_s, h, w_0 and sigma_8:<br/>
```python
   MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}
```
<br/>
You can now compute a non-linear power spectrum at redshift z=0.0 by executing:<br/>

```python
   z = 0.0
   result = e2py.get_pnonlin(MyCosmo, z)
```

The `result` is also a python dictionary with the keys `B` (the boost factor), `P_lin` (the linear power spectrum as computed by class), `P_nonlin` (the non-linear power spectrum being the product of `P_lin` and `B`) and `k` (the vector of k-values given in _h_/Mpc). 

Notice that you can as well pass a list or numpy.ndarray for `z` in which case `result` will contain a nested dictionary: The fields `B`, `P_lin` and `P_nonlin` are now dictionaries themselves with fields `z1`, `z2`, `z3`, etc., corresponding to the boost factor, linear and non-linear power spectrum, respectively, evaluated at the individual redshifts listed in the `z` variable. Here is a short example:
```python
   import e2py
   import numpy as np
   MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}
   z = np.linspace(0.0,1.0,11) # z = [0.0, 0.1, 0.2, ..., 0.9, 1.0]
   result = e2py.get_pnonlin(MyCosmo, z)
```
There will be e.g. the field `result[B][z1]` corresponding to the boost factor at redshift `z=0.0` or the field `result[P_nonlin][z9]` corresponding to the non-linear matter power spectrum evaluated at redshift `z=0.8`. 

## Credits
Credits for this project go to:<br/>
<br/>
__Doug Potter__, University of Zurich, Switzerland (https://bitbucket.org/dpotter/) for constantly advising me in this project <br/><br/>
__Jeppe Mosgaard Dakin__, Aarhus University, Denmark (https://github.com/jmd-dk) for significant contributions to the cython code <br/><br/>
__Rongchuan Zhao__, University of Bonn, Germany (https://astro.uni-bonn.de/m/rzhao) for contributions in the lensing module ee_lens.py and some functions in ee_aux.py.

## License
To be added

## Changes with respect to version 1.0
1. Python wrapper with many useful functions (reducing the risk of falling into a pitfalls)
2. Usage of CMake for the building
