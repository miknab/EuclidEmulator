#!/bin/bash

# ==============
# EXAMPLE SCRIPT
# ============== 

# Executing this script will produce the boost factor for the euclid reference cosmology 
# (as quoted in the Euclid Emulator paper by Knabenhans et al, 2018) at redshifts z1=0.0, z2=0.5, z3=1.0 and z4=2.0 
# and store the result in a data file called "EXAMPLE_EucRefBoost.dat"

ee 0.0219961 0.1431991 0.96 0.67 -1.0 0.83 0.0 0.5 1.0 2.0 > EXAMPLE_EucRefBoost.dat

