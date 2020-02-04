"""
ee_observables.py

EuclidEmulator submodule for actual emulation of cosmological observables.
"""

# This file is part of EuclidEmulator
# Copyright (c) 2018-2020 Mischa Knabenhans
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

import sys as _sys
import numpy as _np
import warnings as _warnings
import EuclidEmulator_BackEnd as _eeb
from e2py._internal import _ee_aux as _aux
from e2py._internal import _ee_background as _bg
from e2py._internal import _ee_cosmoconv as _cc
import e2py.ee_input as _inp
import e2py._ee_lens as _lens

from scipy.integrate import romb as _romb
from scipy.interpolate import CubicSpline as _CubicSpline

try:
    from classy import Class as _Class
except ImportError:
    print "\nClassy could not be found in your system."
    print "Here are some suggestions:\n"
    print "\t -Download the Class from class-code.net and install it"
    print "\t  together with its wrapper classy (type 'make' instead of"
    print "\t  'make class'"
    print "\t -If you know that Class is installed on your system"
    print "\t  and yet classy could not be installed, try re-compiling"
    print "\t  Class with just ''make'' instead of ''make class''"
    print "NOTICE: Even without classy you can still use EuclidEmulator"
    print "        to emulate boost factors. You won't be able to compute"
    print "        full power spectra, though."

def get_boost(emu_pars_dict, redshifts, custom_kvec=None, verbose=True):
    """
    Signature:   get_boost(emu_pars_dict, redshifts [, custom_kvec=None, verbose=True])

    Description: Computes the non-linear boost factor for a cosmology
                 defined in emu_pars_dict (a python dictionary containing
                 the values for the 6 LCDM parameters) at specified
                 redshift stored in a list or numpy.ndarray.
                 Optionally, a list or numpy.ndarray of k modes can be
                 passed to the function via the keyword argument "kvec".
                 Then, by setting verbose=False, it is possible to fully
                 suppress any verbose information about how the code
                 progresses. 

    Input types: python dictionary (with the six cosmological parameters)
                 list or numpy.ndarray (with redshift values)

                 :OPTIONAL:
                 list or numpy.ndarray (with k mode values)
                 boolean (verbosity)

    Output type: python dictionary

    Related:     get_plin, get_pnonlin
    """
    # Check cosmological parameter ranges
    _inp.check_param_range(emu_pars_dict)

    if isinstance(redshifts, (int, float)):
        redshifts = _np.asarray([redshifts])
    else:
        redshifts = _np.asarray(redshifts)

    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    if not isinstance(emu_pars_dict, (dict,)):
        print "The cosmological parameters must be passed as a python \
               dictionary.\n"
        _sys.exit()

    boost_data = _eeb.emu_boost(_np.array([emu_pars_dict['om_b'],
                                           emu_pars_dict['om_m'],
                                           emu_pars_dict['n_s'],
                                           emu_pars_dict['h'],
                                           emu_pars_dict['w_0'],
                                           emu_pars_dict['sigma_8']]),
                                redshifts, verbose)

    kvals = boost_data.k
    k_shape = kvals.shape

    do_extrapolate_above = False
    do_extrapolate_below = False
    if not(custom_kvec is None):
        upper_mask = custom_kvec < max(kvals)
        lower_mask = custom_kvec > min(kvals)
        mask = [u and l for (u,l) in zip(lower_mask, upper_mask)]
        custom_k_within_range = custom_kvec[mask]
        custom_k_below = custom_kvec[[not(l) for l in lower_mask]]
        custom_k_above = custom_kvec[[not(u) for u in upper_mask]]

        if any(custom_kvec > max(kvals)):
            wrn_message = ("EuclidEmulator emulates the non-linear correction in \n"
                           "the interval [6.87215e-3 h/Mpc, 5.52669h/Mpc]. You are \n"
                           "requesting k modes beyond k_max = 5.52669h/Mpc. \n"
                           "Higher k modes constantly extrapolated.")
            if verbose:
                _warnings.warn(wrn_message)
            do_extrapolate_above = True

        if any(custom_kvec < min(kvals)):
            wrn_message = ("EuclidEmulator emulates the non-linear correction in \n"
                           "the interval [6.87215e-3 h/Mpc, 5.52669h/Mpc]. You are \n"
                           "requesting k modes below k_min = 6.87215h/Mpc. \n"
                           "Lower k modes constantly extrapolated.")
            if verbose:
                _warnings.warn(wrn_message)
            do_extrapolate_below = True

    len_kvec = len(kvals)
    len_redshifts = len(redshifts)

    if len_redshifts > 1:
        bvals = {}
        for i in range(len_redshifts):
            tmp = boost_data.boost[i*len_kvec:(i+1)*len_kvec]
            if not(custom_kvec is None):
                bvals['z'+str(i)] = 10.0**_CubicSpline(_np.log10(kvals),
                                                       _np.log10(tmp.reshape(k_shape))
                                                      )(_np.log10(custom_k_within_range))

                #Extrapolate if necessary
                if do_extrapolate_below:
                    # below the k_min of EuclidEmulator, we are in the linear regime where
                    # the boost factor is unity by construction
                    b_extrap = _np.ones_like(custom_k_below)
                    bvals['z'+str(i)]= _np.concatenate((b_extrap, bvals['z'+str(i)]))

                if do_extrapolate_above:
                    # We extrapolate by setting all b(k > k_max) to b(k_max)
                    b_extrap = bvals['z'+str(i)][-1] * _np.ones_like(custom_k_above)
                    bvals['z'+str(i)] = _np.concatenate((bvals['z'+str(i)], b_extrap))

            else:
                bvals['z'+str(i)] = tmp.reshape(k_shape)
    else:
        tmp = boost_data.boost
        if not(custom_kvec is None):
            bvals = 10.0**_CubicSpline(_np.log10(kvals), 
                                       _np.log10(tmp.reshape(k_shape))
                                      )(_np.log10(custom_k_within_range))

            #Extrapolate if necessary
            if do_extrapolate_below:
                # below the k_min of EuclidEmulator, we are in the linear regime where
                # the boost factor is unity by construction
                b_extrap = _np.ones_like(custom_k_below)
                bvals = _np.concatenate((b_extrap,bvals))

            if do_extrapolate_above:
                # We extrapolate by setting all b(k > k_max) to b(k_max)
                b_extrap = bvals[-1] * _np.ones_like(custom_k_above)
                bvals = _np.concatenate((bvals, b_extrap))

        else:
            bvals = tmp.reshape(k_shape)

    if not(custom_kvec is None):       # This could probably be done cleaner!
        kvals = custom_kvec

    return {'k': kvals, 'B': bvals}

def get_pnonlin(emu_pars_dict, redshifts, custom_kvec=None, verbose=True):
    """
    Signature:   get_pnonlin(emu_pars_dict, redshifts [, custom_kvec=None, verbose=True])

    Description: Computes the linear power spectrum and the non-linear boost
                 separately for a given redshift z (or for a list or numpy.ndarray
                 of redshifts), a given cosmology defined in emu_pars_dic (a python
                 dictionary containing the values for the 6 LCDM parameters) and
                 optionally a list or numpy.ndarray of k modes. Then it returns the
                 product of these two which is the non-linear DM-only power spectrum.

    Input types: python dictionary (with the six cosmological parameters)
                 float or iterable (list, numpy.ndarray) (with redshifts)

                 :OPTIONAL:
                 iterable (list, numpy.ndarray) (with k modes)
                 boolean (verbose)

    Output type: python dictionary

    Related:     get_plin, get_boost
    """
    if _Class.__module__ not in _sys.modules:
        print "You have not imported neither classee nor classy.\n \
               Emulating full power spectrum is hence not possible."
        return None

    # Check cosmological parameter ranges
    _inp.check_param_range(emu_pars_dict)

    if isinstance(redshifts, (int, float)):
        redshifts = _np.asarray([redshifts])
    else:
        redshifts = _np.asarray(redshifts)

    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    boost_dict = get_boost(emu_pars_dict, redshifts, custom_kvec, verbose)

    kvec = boost_dict['k']
    Bk = boost_dict['B']

    plin = get_plin(emu_pars_dict, kvec, redshifts)
    plin = plin['P_lin']

    if len(redshifts) == 1:
        pnonlin = plin*Bk

    else:
        pnonlin = {}
        for i, z in enumerate(redshifts):
            pnonlin['z'+str(i)] = plin['z'+str(i)]*Bk['z'+str(i)]

    return {'k': kvec, 'P_nonlin': pnonlin, 'P_lin': plin, 'B': Bk}

def get_plin(emu_pars_dict, custom_kvec, redshifts):
    """
    Signature:   get_plin(emu_pars_dict, custom_kvec, redshifts)

    Description: Computes the linear power spectrum at redshift z for a
                 cosmology defined in EmuParsArr (a numpy array containing
                 the values for the 6 LCDM parameters) (uses classy).

    Input types: python dictionary (with the six cosmological parameters)
                 numpy.ndarray (containing the k modes)
                 numpy.ndarray (containing the redshift values)

    Output type: if len(redshifts)==1, then numpy.ndarray (containing the
                 linear power spectrum values)
                 if len(redshifts)>1, then dict with indices 'z0', 'z1',
                 'z2' etc.

    Related:     get_pnonlin, get_boost
    """
    if _Class.__module__ not in _sys.modules:
        print "You have not imported neither classee nor classy.\n \
               Computing linear power spectrum is hence not possible."
        return None

    # Convert single redshift input argument to array
    if isinstance(redshifts, (int, float)):
        redshifts = _np.asarray([redshifts])
    else:
        redshifts = _np.asarray(redshifts)

    for z in redshifts:
        assert z <= 5.0, "EuclidEmulator allows only redshifts z <= 5.0.\n"

    # Convert single redshift input argument to array
    if isinstance(custom_kvec, (int, float)):
        custom_kvec = _np.asarray([custom_kvec])
    else:
        custom_kvec = _np.asarray(custom_kvec)

    # "Stringify" the input arrays to be understandable for classy.
    z_arr, z_str = _aux.stringify_arr(redshifts)

    # Convert the input dictionary into a Class-compatible dictionary
    class_pars_dict = _inp.emu_to_class(emu_pars_dict)

    # Extend the input Class-compatible dictionary by the additional
    # information requested by classy.
    classy_pars = class_pars_dict
    classy_pars['Omega_Lambda'] = 0.0
    classy_pars['output'] = 'mPk'
    classy_pars['P_k_max_1/Mpc'] = 10.0
    classy_pars['z_pk'] = z_str

    # Create a "Class" instance called "cosmo" and run classy to compute
    # the cosmological quantities.
    cosmo = _Class()
    cosmo.set(classy_pars)
    cosmo.compute()


    # Convert k units: h/Mpc
    h = classy_pars['h']
    k_classy_arr = h*custom_kvec

    # Get shape of k vector
    k_shape = k_classy_arr.shape

    # Get power spectrum at tabulated z and k in units of Mpc^3
    if len(z_arr) == 1:
        z = z_arr
        linpower = _np.array([cosmo.pk(k, z)*h*h*h
                              for k in k_classy_arr]).reshape(k_shape)
    else:
        linpower = {'z'+str(i):
                    _np.array([cosmo.pk(k, z)*h*h*h
                               for k in k_classy_arr]).reshape(k_shape)
                    for i, z in enumerate(z_arr)}

    return {'k': custom_kvec, 'P_lin': linpower}
