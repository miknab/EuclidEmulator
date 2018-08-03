"""
ee_observables.py

EuclidEmulator submodule for actual emulation of cosmological observables.
"""

import sys
import numpy as np
import EuclidEmulator_BackEnd as eeb
import ee_aux as aux
import ee_input as inp
import ee_background as bg
import ee_lens as lens
import ee_cosmoconv as cc

from scipy.integrate import romb
from scipy.interpolate import CubicSpline
from classy import Class

def get_boost(emu_pars_dict, redshifts):
    """
    Signature:   get_boost(emu_pars_dict, redshifts)

    Description: Computes the non-linear boost factor for a cosmology
                 defined in EmuParsArr (a numpy array containing the
                 values for the 6 LCDM parameters) at specified
                 redshift stored in a list or numpy.array.

    Input types: python dictionary (with the six cosmological parameters)
                 list or numpy.array

    Output type: python dictionary

    Related:     get_plin, get_pnonlin
    """
    redshifts = np.asarray(redshifts)
    
    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    if not isinstance(emu_pars_dict, (dict,)):
        print("The cosmological parameters must be passed as a python \
               dictionary.\n")
        sys.exit()

    boost_data = eeb.emu_boost(np.array([emu_pars_dict['om_b'],
                                         emu_pars_dict['om_m'],
                                         emu_pars_dict['n_s'],
                                         emu_pars_dict['h'],
                                         emu_pars_dict['w_0'],
                                         emu_pars_dict['sigma_8']]),
                               redshifts)
    
    kvals = boost_data.k
    k_shape = kvals.shape

    nk = len(kvals)
    nz = len(redshifts)

    if nz>1:
        bvals = {}
        for i in range(nz):
            bvals['z'+str(i)]=boost_data.boost[i*nk:(i+1)*nk].reshape(k_shape) 
    else:
        bvals = boost_data.boost.reshape(k_shape)

    return {'k': kvals, 'B': bvals}

def get_pnonlin(emu_pars_dict, redshifts):
    """
    Signature:   get_pnonlin(emu_pars_dict, redshifts)

    Description: Computes the linear power spectrum and the non-linear boost
                 separately for a given redshift z and cosmology defined by
                 EmuParsArr (a numpy array containing the values for the 6 LCDM
                 parameters) and then returns the product of these two which is
                 the non-linear DM-only power spectrum.

    Input types: python dictionary (with the six cosmological parameters)
                 iterable (list, numpy array)

    Output type: python dictionary

    Related:     get_plin, get_boost
    """
    redshifts = np.asarray(redshifts)

    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"
 
    boost_dict = get_boost(emu_pars_dict, redshifts)

    kvec = boost_dict['k']
    shape = kvec.shape

    Bk = boost_dict['B']

    plin = get_plin(emu_pars_dict, kvec, redshifts)

    if len(redshifts)==1:
        pnonlin = plin*Bk

    else:
        pnonlin = {}
        for i,z in enumerate(redshifts):
            pnonlin['z'+str(i)] = plin['z'+str(i)]*Bk['z'+str(i)]
    
    return {'k': kvec, 'P_nonlin': pnonlin, 'P_lin': plin, 'B': Bk}

def get_plin(emu_pars_dict, k_arr, z_arr):
    """
    Signature:   get_plin(emu_pars_dict, k_arr, z_arr)

    Description: Computes the linear power spectrum at redshift z for a
                 cosmology defined in EmuParsArr (a numpy array containing
                 the values for the 6 LCDM parameters) (uses classy).

    Input types: python dictionary (with the six cosmological parameters)
                 numpy.ndarray (containing the k modes)
                 numpy.ndarray (containing the redshift values)

    Output type: if len(z_arr)==1, then numpy.ndarray (containing the linear power spectrum values)
                 if len(z_arr)>1, then dict with indices 'z0', 'z1', 'z2' etc.

    Related:     get_pnonlin, get_boost
    """
    # Convert single redshift input argument to array
    if isinstance(z_arr, (float, int)):
        z_arr = np.array([z_arr])

    for z in z_arr:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    # "Stringify" the input arrays to be understandable for classy.
    z_arr, z_str = aux.stringify_arr(z_arr)
    k_arr, k_str = aux.stringify_arr(k_arr)

    # Convert the input dictionary into a Class-compatible dictionary
    class_pars_dict = inp.emu_to_class(emu_pars_dict)

    # Extend the input Class-compatible dictionary by the additional
    # information requested by classy.
    classy_pars = class_pars_dict
    classy_pars['Omega_Lambda'] = 0.0
    classy_pars['output'] = 'mPk'
    classy_pars['P_k_max_1/Mpc'] = 10.0
#    classy_pars['k_output_values'] = k_str
    classy_pars['z_pk'] = z_str

    # Create a "Class" instance called "cosmo" and run classy to compute
    # the cosmological quantities.
    cosmo = Class()
    cosmo.set(classy_pars)
    cosmo.compute()

    # Convert k units: h/Mpc
    h = classy_pars['h']
    k_classy_arr = h*k_arr

    # Get shape of k vector
    k_shape = k_classy_arr.shape

    # Get power spectrum at tabulated z and k in units of Mpc^3 
    if len(z_arr)==1:
        z = z_arr
        linpower = np.array([cosmo.pk(k, z)*h*h*h for k in k_classy_arr]).reshape(k_shape)
    else:
        linpower = {'z'+str(i): np.array([cosmo.pk(k, z)*h*h*h for k in k_classy_arr]).reshape(k_shape) for i,z in enumerate(z_arr)}

    return linpower

def get_pconv(emu_pars_dict, sourcedist, prec=12):
    """
    Signature:   get_pconv(emu_pars_dict, sourcedist, prec=12)
        
    Description: Converts a matter power spectrum into a convergence power 
                 spectrum via the Limber equation (conversion is based on
                 equation (29) of the review article by Martin Kilbinger 
                 "Cosmology with cosmic shear observations: a review", July 21,
                 2015,arXiv:1411.0115v2). REMARK: The g in that equation is a
                 typo and should be q which is defined in equation (24). 
                 
                 The input variable emu_pars_dict is a dictionary containing
                 the cosmological parameters and sourcedist is a dictionary of
                 the format {'chi': ..., 'n': ...} containing the a vector of 
                 comoving distances ("chi") together with the number counts of 
                 source galaxies at these distances or redshifts, respectively.
                 
                 REMARK: The 'chi' vector in the sourcedist dictionary should 
                 span the range from 0 to chi(z=5).

    Input type:  emu_pars_dict - dictionary
                 sourcedist - np.ndarray
        
    Ouput type:  dictionary of the form {'l': ..., 'Pconv': ...}
    """
    c = 1 # speed of light
    c_inv = 1./c
    h = emu_pars_dict['h']
    H0 = 100*h
    Om_m = emu_pars_dict['om_m']/(h*h)
    
    nz = int(2**prec + 1)
    nl = int(1e4)
    
    z_vec = np.linspace(0.0, 5.0, nz)
    a_inv_vec = 1.+z_vec
    chi_lim = bg.dist_comov(emu_pars_dict, 0.0, 5.0)
    chi_vec = []
    pnonlin_array = []
    
    l_vec = np.linspace(1,1e3, nl)
    
    for z in z_vec:
        chi_vec.append(bg.dist_comov(emu_pars_dict, 0.0, z))
        P = get_pnonlin(emu_pars_dict, z) # call EuclidEmulator
        
        print len(P['k']), len(P['P_nonlin'])
        
        f = CubicSpline(np.log10(P['k']), np.log10(P['P_nonlin']))
        k = cc.l_to_k(emu_pars_dict, l_vec, z)
        pnonlin_array.append(10** f(np.log10(k)))
    
    chi_vec = np.asarray(chi_vec)
    pnonlin_array = np.asarray(pnonlin_array)
        
    prefac = 2.25*Om_m*Om_m * H0 * H0 * H0 * H0 * c_inv * c_inv * c_inv * c_inv
    
    # compose the integrand
    assert (len(z_vec) == len(chi_vec))
    assert (len(chi_vec) == len(a_inv_vec))
    assert all([len(chi_vec) == len(pnonlin_array[i]) 
                for i in range(len(pnonlin_array))])
    
    q_vec = lens.lens_efficiency(sourcedist, chi_vec, chi_lim)
 
    integrand = q_vec * q_vec * a_inv_vec * a_inv_vec * pnonlin_array

    assert(integrand.shape == pnonlin_array.shape)
    
    # perform integral and return
    Pconv = {'l': l_vec, 
             'Pconv': prefac*romb(integrand, chi_vec[1]-chi_vec[0])}
    
    return Pconv
