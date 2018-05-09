# EuclidEmulator
This repository contains the main source code of the EuclidEmulator, a fast and accurate tool to estimate the non-linear correction to the matter power spectrum.

Authors:   M. Knabenhans & J. Stadel<br/>
Date:      May 2018<br/>
Reference: Knabenhans et al., 2018, arXiv pre-print (submitted)<br/>

## Description
EuclidEmulator is a tool for rapid and accurate estimation of the
non-linear correction ("boost") to the dark matter power spectrum.
It is based on a spectral decomposition technique (called "sparse
polynomial chaos expansion") and estimates the boost at an accuracy
level of a fraction of a percent in the entire allowed redshift and
wavenumber range.

An estimate of a full non-linear power spectrum is obtained by
multiplying the emulated boost with a linear power spectrum computed
by a Boltzmann code (e.g. CAMB or CLASS). As a result, the linear
regime of the non-linear power spectrum is quasi perfect and effects
not captured by Newtonian gravity of dark matter (such as super-horizon
damping or baryonic clustering) are accounted for to first order on
all scales.

For more details we refer to publication.

## Parameter ranges
For the emulation of the boost, the cosmological parameters have to be
within the following ranges:

0.0217 <= om_b    <= 0.0233<br/>
0.1326 <= om_m    <= 0.1526<br/>
0.9345 <= n_s     <= 0.9965<br/>
0.6251 <= h       <= 0.7211<br/>
-1.250 <= w_0     <= -0.750<br/>
0.7684 <= sigma_8 <= 0.8614<br/>

The redshift has to be 0 <= z <= 5.

Input values outside this range will produce an error.

## Code structure
### Main emulator code
The emulator source code consists of four files:

`ee.dat` (!!! NEVER EVER CHANGE THIS FILE !!!)<br/>
`ee.c`<br/>
`cosmo.c`<br/>
`cosmo.h`<br/>

In ee.dat, all the building blocks of the emulator are stored. This file
forms the core of the emulator code. The numbers in ee.dat are the results
of the study described in Knabenhans et al. (2018). It includes the mean
boost factor, the principal components, the normalized Legendre polynomials
and the corresponding coefficients. MODIFYING THIS FILE WILL RESULT IN WRONG
EMULATION RESULTS.

The code in ee.c is just the front end that reads in, assembles and evaluates
the functions and coefficients stored in ee.dat. From this file, the executable
"ee" is built via the Makefile.

The code stored in cosmo.c and cosmo.h, respectively, is used for the inter-
polation in redshift space (in order to allow the user to ask for any real-
valued redshift z <= 5).

### Auxiliary scripts and files
In addition, the following files are in the EuclidEmulator repository:

`EuclidRef_Class.ini`<br/>
`GetPnonlin.py`<br/>
`example.sh`<br/>

`EuclidRef_Class.ini` is a CLASS parameter file specifying the relevant parameters for the Euclid reference cosmology. REMARK: \Omega_{rad} is uniquely determined by the CMB temperature. Since \Omega_{rad} was fixed for all cosmologies in the experimental design (cf. reference paper by Knabenhans et al., 2018; submitted), it must not be changed. Notice further, that the parameter A_s has to be chosen in accordance to the Emulator input parameter \sigma_8. 

The python script `GetPnonlin.py` 

## User Guide
1. Prerequisites:<br/>
   * GNU Scientific Library (GSL)
   * Python2.7 (for post-processing only)
   * CLASS (for post-processing only)
   
2. Building the code:<br/>
   Type "make". An executable called "ee" will be created.

   Make sure not to move the files cosmo.c, cosmo.h and ee.dat.
   They have to be in the same directory as EuclidEmulator.c and
   the executable ee.

3. Usage:<br/>
   ./ee <om_b> <om_m> <n_s> <h> <w_0> <sigma_8> <z>

   This will print the resulting boost factor to standard output. To store
   it in a file, just use output redirection, i.e. append " > BoostFile.dat"
   to the command above.

4. Computation of the full non-linear power spectrum:<br/>
   1. Produce a linear power spectrum with a Boltzmann code.
   2. Interpolate the linear power spectrum (we suggest a cubic spline interpolation in log space.
   3. Evaluate the interpolated liner power spectrum at the k modes of the emulated boost.
   4. Multiply the boost with the resulting linear power spectrum.

The script `example.sh` includes commands to perform step iv (assuming you have CLASS installed already and your CLASS file explanatory.ini has never been changed). Please change the path to the directory where your CLASS executable is located. Executing `example.sh` will sequentially compute the boost spectrum using EuclidEmulator, compute the corresponding linear power spectrum with CLASS, and finally call `GetPnonlin.py` to compute the non-linear power spectrum.
