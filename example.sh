#!/bin/bash

# ==============
# EXAMPLE SCRIPT
# ============== 

# Executing this script will produce the boost factor for the euclid reference cosmology 
# (as quoted in the Euclid Emulator paper by Knabenhans et al, 2018) at redshift z=0.5 
# and store the result in a data file called "EXAMPLE_EucRefBoost.z0.5.dat"

./ee 0.0219961 0.1431991 0.96 0.67 -1.0 0.83 0.5 > EXAMPLE_EucRefBoost.z0.5.dat

/home/ics/mischak/MyApplications/CLASS/class_public-2.6.3/class EucRef_Class.ini

python GetPnonlin.py
