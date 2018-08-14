# EuclidEmulator (version 1.1)
This repository contains the main source code of the EuclidEmulator, a fast and accurate tool to estimate the non-linear correction to the matter power spectrum.

Authors:   M. Knabenhans & J. Stadel<br/>
Date:      May 2018<br/>
Reference: Knabenhans et al., 2018, arXiv pre-print (submitted)<br/>

If you use EuclidEmulator in any way (for a publication or otherwise), please cite this paper.

STAY TUNED:
1) We are working on a new feature for the python wrapper allowing to easily and quickly calculate convergence power spectra.
2) We plan a fully revised version of the emulator that includes neutrino physics.

## Quick start
If you have not done so already, either download this repository or clone it to your local host (under Linux you can get a gzipped tar-ball via `wget https://github.com/miknab/EuclidEmulator/archive/master.tar.gz`). In case you have to do so, extract the archive and enter it.

Depending on whether you want to build the python wrapper (recommended) or the command line interface (CLI) only, enter one or the other directory (called `wrapper` or `CLI`). Create a new build directory (in bash e.g. via `mkdir build`), enter this newly created directory and type <br/><br/>
`cmake ..`  <br/><br/>
(if you have root access) or <br/><br/>
`cmake -DCMAKE_INSTALL_PREFIX=path/to/installation/directory ..`  <br/><br/>
(without root access). Then type `make` in order to build the code and `make install` to install it.

If you don't have root acces on the machine you want to build this software on, we recommend to install this software inside a virtual environment (see documentation/wiki for more info).

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

For more details we refer to publication and the documentation/wiki.

## Changes with respect to version 1.0
1. Python wrapper with many useful functions (reducing the risk of falling into a pitfall)
2. Usage of CMake for the building

## User guide

=============================== I M P O R T A N T ===============================<br/>
!!! PLEASE READ THE SECTION ABOUT PITFALLS BELOW IN ORDER TO ENSURE CORRECT RESULTS !!!
============================================================================<br/>

1. Prerequisites:<br/>
   * GNU Scientific Library (GSL)
   * Python2.7 together with numpy, scipy and matplotlib (for post-processing only)
   * CLASS (for post-processing only)
   
2. Building the code:<br/>
   Type "make". An executable called "ee" will be created.

   Make sure not to move the files cosmo.c, cosmo.h and ee.dat.
   They have to be in the same directory as EuclidEmulator.c and
   the executable ee.

3. Usage:<br/>
   ./ee <&#969;<sub>b</sub>> <&#969;<sub>m</sub>> <n<sub>s</sub>> \<h\> <w<sub>0</sub>> <&#963;<sub>8</sub>> \<z\>

   This will print the resulting boost factor to standard output. To store
   it in a file, just use output redirection, i.e. append " > BoostFile.dat"
   to the command above.

4. Computation of the full non-linear power spectrum:<br/>
   1. Produce a linear power spectrum with a Boltzmann code.
   2. Interpolate the linear power spectrum (we suggest a cubic spline interpolation in log space.
   3. Evaluate the interpolated liner power spectrum at the k modes of the emulated boost.
   4. Multiply the boost with the resulting linear power spectrum.

The script `example.sh` includes commands to perform step iv (assuming you have CLASS installed already and your CLASS file explanatory.ini has never been changed). Please change the path to the directory where your CLASS executable is located. Executing `example.sh` will sequentially compute the boost spectrum using EuclidEmulator, compute the corresponding linear power spectrum with CLASS, and finally call `GetPnonlin.py` to compute the non-linear power spectrum.

## Pitfalls

When you want to emulate a full non-linear power spectrum, you really have to make sure that you specify the exact same cosmology for EuclidEmulator to produce the boost factor and for the Boltzmann code you use to predict the linear power spectrum. Often, the parametrization of cosmologies used in the Boltzmann solvers is different than the one used by EuclidEmulator. Make sure the different parametrizations define the exact same cosmology!

Known differences are:

1.  CAMB and CLASS use om_cdm (the cold dark matter density &#969;<sub>cdm</sub>) instead of om_m (the total matter density
    &#969;<sub>m</sub>). Make sure that the following relation is satisfied: 
    <div align="center">&#969<sub>b</sub> + &#969<sub>cdm</sub> = &#969<sub>m</sub> </div>
    If you want to compute the non-linear power spectrum using a EuclidEmulated-boost for a real data set (like e.g. 
    Planck2015), you may find that the reported values doe not obey the above relation. This is most likely due to the fact
    that there is a small contribution due to neutrinos which is taken into account in that data set but not so in this 
    version of our emulator. Here's what you need to do:</br>

    1. To produce the boost factor with EuclidEmulator use &#969;<sub>b</sub> and &#969;<sub>m</sub> as reported in the data set.
    2. To produce the corresponding linear power spectrum use &#969;<sub>b</sub> and set the &#969;<sub>cdm</sub> parameter equal to &#969;<sub>m</sub>-&#969;<sub>b</sub>. Doing so you add the neutrino component to the CDM contribution which is the best that can be done with the current version of EuclidEmulator. Stay tuned as version 2 will allow for neutrinos to be taken into accound.
 
2.  CAMB and CLASS do usually not accept sigma_8 as a parameter for normalization of the power spectrum but rather use A_s. In
    order to convert these two parameters into each other in the context of using EuclidEmulator, you have to use the same
    conversion as is used in the EuclidEmulator code. Convert the parameters using the following proportionality:<br/>
    <div align="center"> A<sub>s</sub>/(2.215 * 10^(-9)) = (&#963;<sub>8</sub>/0.8496)^2
