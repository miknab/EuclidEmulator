"""
ee_background.py

EuclidEmulator submodule for computation of cosmological background quantities.

REMARK:      The geometry of the Universe is fixed to be flat (i.e.
             Omega_curvature = 1) and the radiation energy density
             is set to Om_rad = 4.183709411969527e-5/(h*h). These
             values were assumed in the construction process of 
             EuclidEmulator and hence must be used whenever one is
             working with it.
"""

# This file is part of EuclidEmulator 
# Copyright (c) 2018 Mischa Knabenhans
#
# EuclidEmulator is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# EuclidEmulator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy.integrate import romb
import _ee_cosmoconv as cc

# UNIVERSAL CONSTANTS (Source: PDG booklet 2018)
SPEED_OF_LIGHT_IN_KILOMETERS_PER_SECOND = 299792.458
NEWTONS_CONSTANT = 6.6740908 * 1e-11 # in m^3kg^(-1)s^(-2)
MEGAPARSEC = 3.08567758149 * 1e22 # meters
M_SOLAR_IN_KG = 1.98848 *1e30

def dist_comov(emu_pars_dict, z1, z2, prec=12):
    """
    Signature:   dist_comov(emu_pars_dict, z1, z2, prec=12)

    Description: Computes the comoving distance (in units of Mpc/h) between 
                 objects at redshifts z1 and z2 for a cosmology specified by
                 the parameters Om_m (matter density parameter),
                 Om_rad (radiation density parameter) and Om_DE (dark
                 energy density parameter), the Hubble parameter H0,
                 and the dark energy equation of state parameters w0
                 and wa.

    Input type:  python dictionary (containing the 6 LCDM parameters)
                 floats or np.ndarrays 

    Output type: float or np.ndarray

    REMARK 1:    If the redshifts z1 and z2 are passed as vectors (1-dimen-
                 sional np.ndarrays) then the comoving distances will be com-
                 puted between the pairs of redshift (z1_1, z2_1), (z1_2,z2_2),
                 (z1_3, z2_3) etc. Hence, if the vectors have length n, the 
                 resulting vector of d_comov will also be of length n (and not
                 n(n-1)). We do NOT compute the d_comov for all possible red-
                 shift combinations.

    REMARK 2:    The geometry of the Universe is fixed to be flat (i.e.
                 Omega_curvature = 1) and the radiation energy density
                 is set to Om_rad = 4.183709411969527e-5/(h*h). These
                 values were assumed in the construction process of 
                 EuclidEmulator and hence must be used whenever one is
                 working with it.
    """
       
    if isinstance(z1,(float,int)) and isinstance(z2,(float,int)):
        z1 = np.array([z1])
        z2 = np.array([z2])

    a1_vec = cc.z_to_a(z1)
    a2_vec = cc.z_to_a(z2)

    d_comov = []
    for a1,a2 in zip(a1_vec,a2_vec):
        if a1 > a2:
            a1, a2 = a2, a1 # swap values
        avec = np.linspace(a1, a2, 2**prec+1)
        delta_a = avec[1]-avec[0]

        H = cc.a_to_hubble(emu_pars_dict,avec)

        d_comov.append(romb(1./(avec*avec*H), dx=delta_a))

    chi = np.array(d_comov) # here the comoving distances have units of Mpc/(km/s)
    
    # return result in units of Mpc/h
    return SPEED_OF_LIGHT_IN_KILOMETERS_PER_SECOND * chi * emu_pars_dict['h']
