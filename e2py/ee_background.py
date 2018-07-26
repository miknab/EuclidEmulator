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

import numpy as np
from scipy.integrate import romb
import ee_cosmoconv as cc

def dist_comov(emu_pars_dict, z1, z2, prec=12):
    """
    Signature:   dist_comov(emu_pars_dict, z1, z2, prec=12)

    Description: Computes the comoving distance between objects at
                 redshifts z1 and z2 for a cosmology specified by
                 the parameters Om_m (matter density parameter),
                 Om_rad (radiation density parameter) and Om_DE (dark
                 energy density parameter), the Hubble parameter H0,
                 and the dark energy equation of state parameters w0
                 and wa.

    Input type:  python dictionary (containing the 6 LCDM parameters)
                 float (redshift) 

    Output type: float

    REMARK:      The geometry of the Universe is fixed to be flat (i.e.
                 Omega_curvature = 1) and the radiation energy density
                 is set to Om_rad = 4.183709411969527e-5/(h*h). These
                 values were assumed in the construction process of 
                 EuclidEmulator and hence must be used whenever one is
                 working with it.
    """

    a1 = cc.z_to_a(z1)
    a2 = cc.z_to_a(z2)

    avec = np.linspace(a1, a2, 2**prec+1)
    delta_a = avec[1]-avec[0]

    H = cc.a_to_hubble(emu_pars_dict,avec)

    d_comov = abs(romb(H, dx=delta_a))

    return d_comov
