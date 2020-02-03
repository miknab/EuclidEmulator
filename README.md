# EuclidEmulator (version 1.2)
This repository contains the main source code of the EuclidEmulator, a fast and accurate tool to estimate the non-linear correction to the matter power spectrum. EuclidEmulator is roughly seven orders of magnitude faster than an N-body simulation that yields results at the same level of accuracy.

Authors:   M. Knabenhans & J. Stadel<br/>
Date of last update:      September 2018<br/>
Reference: Knabenhans et al. (2019), MNRAS, 484, 5509-5529 (https://arxiv.org/abs/1809.04695)<br/>

If you use EuclidEmulator in any way (for a publication or otherwise), please cite this paper.

STAY TUNED:
1) We are working on a new feature for the python wrapper allowing to easily and quickly calculate convergence power spectra.
2) We plan a fully revised version of the emulator that includes neutrino physics.

## Contact information
If you have any questions and/or remarks related to this work, please do not hesitate to send me an email (mischakATphysik.uzh.ch)

## Quick start
### Prerequisites
In any case you need:
 * GNU Scientific Library (GSL; see https://www.gnu.org/software/gsl/)
 * CMake (see https://cmake.org)

In order to use the python wrapper, you need in addition:
 * Python together with cython, numpy, scipy and pandas (only for python wrapper)
 * CLASS code (see http://class-code.net) version 2.6.3 or higher together with the python wrapper classy (only for python wrapper)
 
### Get the code
If you have not done so already, either download this repository or clone it to your local host (under Linux you can get a gzipped tar-ball via
```
   wget https://github.com/miknab/EuclidEmulator/archive/master.tar.gz
```

### Building and installation
EuclidEmulator comes with a command line interface (CLI) and a full-fledged python package. You can choose whether you want to install only one of these tools or both. We recommend to install the python package as it's usage is less error-prone than that of the CLI and it offers much more utilities. In the following we will restrict our description on the python package (if you are interested in the CLI, please consult the wiki).
Now, inside the `EuclidEmulator` directory type
```bash
   cd wrapper2
```
if you want to install the python2 version of the EuclidEmulator wrapper, whereas if you want the python3 version you enter the `wrapper3` directory like
```bash
   cd wrapper3
```
All the subsequent steps are identical for both python versions. Create a new build folder, enter it and run CMake, i.e.: 
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
Now you are ready to use the software. In a python script or in a jupyter notebook (with a kernel of the corresponding python version) you can import the package via <br/>
```python
   import e2py
```
Next, define a cosmology as a python dictionary with the keys `om_b`, `om_m`, `n_s`, `h`, `w_0` and `sigma_8`:<br/>
```python
   MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}
```
<br/>
You can now compute a non-linear power spectrum at redshift z=0.0 by executing:<br/>

```python
   z = 0.0
   result = e2py.get_pnonlin(MyCosmo, z)
```

The `result` is also a python dictionary with the keys `B` (the boost factor, units: dimensionless), `P_lin` (the linear power spectrum as computed by class, units: (Mpc/_h_)<sup>3</sup>), `P_nonlin` (the non-linear power spectrum being the product of `P_lin` and `B`Mpc/_h_)<sup>3</sup>) and `k` (the vector of k-values given in _h_/Mpc). 
The `result` is also a python dictionary with the keys `B` (the boost factor, units: dimensionless), `P_lin` (the linear power spectrum as computed by class, units: (Mpc/_h_)<sup>3</sup>), `P_nonlin` (the non-linear power spectrum, also in units of (Mpc/_h_)<sup>3</sup>, being the product of `P_lin` and `B`) and `k` (the vector of k-values given in _h_/Mpc). 

Notice that you can as well pass a list or numpy.ndarray for `z` in which case `result` will contain a nested dictionary: The fields `B`, `P_lin` and `P_nonlin` are now dictionaries themselves with fields `z1`, `z2`, `z3`, etc., corresponding to the boost factor, linear and non-linear power spectrum, respectively, evaluated at the individual redshifts listed in the `z` variable. Here is a short example:
```python
   import e2py
   import numpy as np
   MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}
   z = np.linspace(0.0,1.0,11) # z = [0.0, 0.1, 0.2, ..., 0.9, 1.0]
   result = e2py.get_pnonlin(MyCosmo, z)
```
There will be e.g. the field `result[B][z1]` corresponding to the boost factor at redshift `z=0.0` or the field `result[P_nonlin][z9]` corresponding to the non-linear matter power spectrum evaluated at redshift `z=0.8`. 

*NOTICE:* The number of redshifts at which you ask EuclidEmulator to produce non-linear power spectra is limited. Usually this upper limit is 100 (sometimes already a bit less). This is a problem of CLASS (there is no such limit for boost factors whose computation does not rely on CLASS).

## Credits
Credits for this project go to:<br/>
<br/>
__Doug Potter__, University of Zurich, Switzerland for constantly advising me in this project; https://bitbucket.org/dpotter/ <br/><br/>
__Jeppe Mosgaard Dakin__, Aarhus University, Denmark for significant contributions to the cython code; https://github.com/jmd-dk <br/><br/>
__Rongchuan Zhao__, University of Bonn, Germany for contributions in the lensing module ee_lens.py and some functions in ee_aux.py; https://astro.uni-bonn.de/m/rzhao.

## License
EuclidEmulator is free software, distributed under the GNU General Public License. This implies that you may freely distribute and copy the software. You may also modify it as you wish, and distribute these modified versions. You must always indicate prominently any changes you made in the original code and leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. 

## Changes with respect to version 1.0
1. Python wrapper with many useful functions (reducing the risk of falling into a pitfalls)
2. Usage of CMake for the building
