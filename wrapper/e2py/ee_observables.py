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

try:
    from classee import Class
except ImportError:
    try:
        from classy import Class
    except ImportError:
        print("Classy could not be found in your system. Here are some suggestions:\n")
        print("\t -Download the patched version of Class and its wrapper classee (see https://github.com/miknab/ClassPatch)")
        print("\t -If you know that Class is installed on your system and yet classy could not be installed, try re-compiling Class with just ''make'' instead of ''make class''")
        print("\t -If both of the two previous options are not applicable to your situation, please download and build the original class code (see class-code.net)")
        print("")
        print("NOTICE: Even without classy you can still use EuclidEmulator to emulate boost factors. You won't be able to compute full power spectra, though.")

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
    if Class.__module__ not in sys.modules:
        print("You have not imported neither classee nor classy. Emulating full power spectrum is hence not possible.")
        return None

    if isinstance(redshifts,(int,float)):
        redshifts = np.asarray([redshifts])
    else:
        redshifts = np.asarray(redshifts)   
 
    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"
 
    boost_dict = get_boost(emu_pars_dict, redshifts)

    kvec = boost_dict['k']
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
    if Class.__module__ not in sys.modules:
        print("You have not imported neither classee nor classy. Computing linear power spectrum is hence not possible.")
        return None

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

def sgaldist(alpha=2.0, beta=1.5, z_mean=0.9):
    return lens.GalaxyRedshiftDist(alpha, beta, z_mean)(z_mean)

def get_pconv(emu_pars_dict, sourcedist_func, prec=7):
    """
    Signature:   get_pconv(emu_pars_dict, sourcedist_func, prec=12)
        
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
                 
                 The precision parameter prec is an integer which defines at
                 how many redshifts the matter power spectrum is emulated for 
                 integration in the Limber equation. The relation between the 
                 number of redshifts nz and the parameter prec is given by 
                 
                                         nz = 2^prec + 1
                 
                 REMARK: The 'chi' vector in the sourcedist dictionary should 
                 span the range from 0 to chi(z=5).

    Input type:  emu_pars_dict - dictionary
                 sourcedist - function object
                 prec - int
                 
    Ouput type:  dictionary of the form {'l': ..., 'Cl': ...}
    """
    if Class.__module__ not in sys.modules:
        print("You have not imported neither classee nor classy. Emulating convergence power spectrum is hence not possible.")
        return None


    c = bg.SPEED_OF_LIGHT_IN_KILOMETERS_PER_SECOND # speed of light
    c_inv = 1./c
    h = emu_pars_dict['h']
    H0 = 100*h
    Om_m = emu_pars_dict['om_m']/(h*h)
    
    nz = int(2**prec + 1)
    nl = int(1e4)
    
    z_vec = np.logspace(np.log10(5e-2), np.log10(4.999999), nz)
    a_inv_vec = 1.+z_vec
    chi_lim = bg.dist_comov(emu_pars_dict, 1e-12, 5.0) # in Mpc/h
    chi_vec = []
    pnonlin_array = []
    
    l_vec = np.linspace(1e1,2e3, nl)

    P = get_pnonlin(emu_pars_dict, z_vec) # call EuclidEmulator
    chi_vec = bg.dist_comov(emu_pars_dict, np.zeros_like(z_vec), z_vec)

    for i,z in enumerate(z_vec):
        f = CubicSpline(np.log10(P['k']), np.log10(P['P_nonlin']['z'+str(i)]))
        k = cc.l_to_k(emu_pars_dict, l_vec, z) # needs to be called inside loop
                                               # because different z-values
                                               # lead to different results

        # evaluate the interpolating function of Pnl for all k in the range
        # allowed by EuclidEmulator (this range is given by P['k'])
        pmatternl = [10.0**f(np.log10(kk)) for kk in k if (kk >= P['k'].min() and kk <= P['k'].max())]
        # for k values below the lower bound or above the upper bound of 
        # this k range, set the contributions to 0.0
        p_toosmall = [0.0 for kk in k if kk < P['k'].min()]
        p_toobig = [0.0 for kk in k if kk > P['k'].max()]

        pnonlin_array.append(np.array(p_toosmall + pmatternl + p_toobig))

    # get the non-linear matter power spectrum in units of [(Mpc/h)^3]
    pnonlin_array = np.asarray(pnonlin_array).transpose()
        
    # compute prefactor of limber equation integral --> in units of [Mpc^4]
    prefac = 2.25*Om_m*Om_m * H0 * H0 * H0 * H0 * c_inv * c_inv * c_inv * c_inv
    
    # convert prefactor units to [(Mpc/h)^4]
    prefac = prefac * h * h * h * h

    # compose the integrand
    assert (len(z_vec) == len(chi_vec))
    assert (len(chi_vec) == len(a_inv_vec))
    assert all([len(chi_vec) == len(pnonlin_array[i]) 
                for i in range(len(pnonlin_array))])
    
    source_dict = {'chi': chi_vec, 'n': sourcedist_func(z_vec)}
    
    q_vec = lens.lens_efficiency(source_dict, chi_vec, chi_lim)
 
    integrand = np.array([q_vec * q_vec * a_inv_vec * a_inv_vec * 
                          pnonlin_array[i] for i in range(nl)])

    assert(integrand.shape == pnonlin_array.shape)
    
    # perform integral (bear in mind that only the product of the 
    # integral and the prefac is truely dimensionless) and return
    Pconv = {'l': l_vec, 
             'Cl': prefac*romb(integrand, chi_vec[1]-chi_vec[0])}
    
    return Pconv
