"""
ee_observables.py

EuclidEmulator submodule for actual emulation of cosmological observables.
"""

import sys
import numpy as np
import EuclidEmulator_BackEnd as eeb
import ee_aux as aux
import ee_input as inp

from classy import Class

def get_boost(emu_pars_dict, z):
    """
    Signature:   get_boost(emu_pars_dict, z)

    Description: Computes the non-linear boost factor for a cosmology
                 defined in EmuParsArr (a numpy array containing the
                 values for the 6 LCDM parameters) at a specified
                 redshift z.

    Input types: python dictionary (with the six cosmological parameters)
                 float

    Output type: python dictionary

    Related:     get_plin, get_pnonlin
    """
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
                                         emu_pars_dict['sigma_8']]), z)

    return {'k': boost_data.k, 'B': boost_data.boost}

def get_pnonlin(emu_pars_dict, z):
    """
    Signature:   get_pnonlin(emu_pars_dict, z)

    Description: Computes the linear power spectrum and the non-linear boost
                 separately for a given redshift z and cosmology defined by
                 EmuParsArr (a numpy array containing the values for the 6 LCDM
                 parameters) and then returns the product of these two which is
                 the non-linear DM-only power spectrum.

    Input types: python dictionary (with the six cosmological parameters)
                 float

    Output type: python dictionary

    Related:     get_plin, get_boost
    """

    assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    boost_dict = get_boost(emu_pars_dict, z)

    kvec = boost_dict['k']
    Bk = boost_dict['B']

    plin = get_plin(emu_pars_dict, kvec, z)

    pnonlin = plin*Bk

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

    Output type: numpy.ndarray (containing the linear power spectrum values)

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
    # Get power spectrum at tabulated z and k in units of Mpc^3
    linpower = np.array([[cosmo.pk(k, z) for k in k_classy_arr] for z in z_arr])

    # Convert P(k) units: Mpc^3 --> Mpc^3/h^3
    #return linpower*h*h*h
    return linpower*(h*h*h)

def get_pconv(k, pdelta, z):
    """
    Converts a matter power spectrum into a convergence power spectrum via
    the Limber equation.
    """
