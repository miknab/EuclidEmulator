"""
ee_lens.py

EuclidEmulator submodule for computation of cosmological lensing quantities.

REMARK:      The geometry of the Universe is fixed to be flat (i.e.
             Omega_curvature = 1) and the radiation energy density
             is set to Om_rad = 4.183709411969527e-5/(h*h). These
             values were assumed in the construction process of 
             EuclidEmulator and hence must be used whenever one is
             working with it.
"""

import numpy as np
from scipy.integrate import romb
from scipy.interpolate import CubicSpline
import ee_cosmoconv as cc

def lens_efficiency(sourcedist, dcomov, dcomov_lim):
    """
    Signature:    lens_efficiency(sourcedist, dcomov, dcomov_lim)
    
    Description:  Computes the lens efficiency function q (see e.g. equation 
                  (24)the review article by Martin Kilbinger "Cosmology with
                  cosmic shear observations: a review", July 21, 2015, 
                  arXiv:1411.0115v2), given a source distribution function n and
                  two comoving distances.
                 
    Input types:  sourcedist - np.ndarray
                  dcomov - float (or int) or np.ndarray
                  dcomov_lim - float
    
    Output types: float (or np.ndarray if type(dcomov)=np.ndarray)
    """    
    
    # interpolate the source distribution function
    nfunc = CubicSpline(np.log10(sourcedist['chi']), np.log10(sourcedist['n']))
    
    result = []
    
    if isinstance(dcomov, np.ndarray):

        for d in dcomov:
            chi = np.linspace(d, dcomov_lim, 1000)
            n = 10**nfunc(np.log10(chi))
            integrand = n * (1-d/chi)
            result.append(romb(integrand, chi[1]-chi[0]))
        
        return np.asarray(result)
    
    elif isinstance(dcomov, float) or isinstance(dcomov,int):
        chi = np.linspace(dcomov, dcomov_lim, 1000)
        n = 10**nfunc(np.log10(chi))
        integrand = n * (1-dcomov/chi)
        
        return romb(integrand, chi[1]-chi[0])
    
    else:
        raise(TypeError, "The second argument 'dcomov' must be either a float,\
                          an integer or a np.ndarray.\n")
        