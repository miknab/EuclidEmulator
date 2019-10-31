"""
ee_observables.py

EuclidEmulator submodule for actual emulation of cosmological observables.
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
    print("\nClassy could not be found in your system.")
    print("Here are some suggestions:\n")
    print("\t -Download the Class from class-code.net and install it")
    print("\t  together with its wrapper classy (type 'make' instead of")
    print("\t  'make class'")
    print("\t -If you know that Class is installed on your system")
    print("\t  and yet classy could not be installed, try re-compiling")
    print("\t  Class with just ''make'' instead of ''make class''")
    print("NOTICE: Even without classy you can still use EuclidEmulator")
    print("        to emulate boost factors. You won't be able to compute")
    print("        full power spectra, though.")

def get_boost(emu_pars_dict, redshifts, custom_kvec=None, verbose=True):
    """
    Signature:   get_boost(emu_pars_dict, redshifts [, custom_kvec=None, verbose=True])

    Description: Computes the non-linear boost factor for a cosmology
                 defined in emu_pars_dict (a python dictionary containing
                 the values for the 6 LCDM parameters) at specified
                 redshift stored in a list or numpy.ndarray. Optionally, 
                 a list or numpy.ndarray of k modes can be passed to the
                 function via the keyword argument "kvec".

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

    if not isinstance(emu_pars_dict, dict):
        print("The cosmological parameters must be passed as a python \
               dictionary.\n")
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
            _warnings.warn(wrn_message)
            do_extrapolate_above = True

        if any(custom_kvec < min(kvals)):
            wrn_message = ("EuclidEmulator emulates the non-linear correction in \n"
                           "the interval [6.87215e-3 h/Mpc, 5.52669h/Mpc]. You are \n" 
                           "requesting k modes below k_min = 6.87215h/Mpc. \n"
                           "Lower k modes constantly extrapolated.")
            _warnings.warn(wrn_message)
            do_extrapolate_below = True

    len_kvals = len(kvals)
    len_redshifts = len(redshifts)

    if len_redshifts > 1:
        bvals = {}
        for i in range(len_redshifts):
            tmp = boost_data.boost[i*len_kvals:(i+1)*len_kvals]
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
                 list or numpy.ndarray (with k mode values)
                 boolean (verbosity)

    Output type: python dictionary

    Related:     get_plin, get_boost
    """
    if _Class.__module__ not in _sys.modules:
        print("You have not imported neither classee nor classy.\n \
               Emulating full power spectrum is hence not possible.")
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
        print("You have not imported neither classee nor classy.\n \
               Computing linear power spectrum is hence not possible.")
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

def sgaldist(alpha=2.0, beta=1.5, z_mean=0.9):
    """
    Signature:      sgaldist(alpha=2.0, beta=1.5, z_mean=0.9)

    Description:    This function takes three (optional) keyword arguments
                    (or any subset thereof) and returns a source galaxy
                    distribution according to the standard distribution function
                    given by

                       n(z;alpha,beta,z_mean) = (z/z_0)^alpha exp[-(z/z_0)^beta]

                    where z_0 = z_mean/1.412. Notice that the default values of
                    alpha, beta and z_mean are as listed in 'Signature'.

    Input:          the input parameters are the parameters of the distribution
                    function n(z) mentioned above.

    Input types:    type(alpha)=float,
                    type(beta)=float,
                    type(z_mean)=float

    Output:         This function returns an object that is similar to a
                    function object in the sense that it takes arguments.
                    To be precise, the returned object takes a numpy.ndarray of
                    redshifts as arguments and computes the correspondingsource
                    galaxy distribution in redshift based on the passed
                    parameters.

    Output types:   type(output)=e2py._internal._ee_aux.Function
                    (EuclidEmulator internal type)

    """
    return _lens.GalaxyRedshiftDist(alpha, beta, z_mean)(z_mean)

if False:
    def get_pconv(emu_pars_dict, sourcedist_func, prec=7):
        """
        Signature:   get_pconv(emu_pars_dict, sourcedist_func, prec=12)

        Description: Converts a matter power spectrum into a convergence power
                     spectrum via the Limber equation (conversion is based on
                     equation (29) of the review article by Martin Kilbinger
                     "Cosmology with cosmic shear observations: a review",
                     July 21, 2015,arXiv:1411.0115v2).
                     REMARK: The g in that equation is a typo and should be q
                     which is defined in equation (24).

                     The input variable emu_pars_dict is a dictionary contai-
                     ning the cosmological parameters and sourcedist is a
                     dictionary of the format {'chi': ..., 'n': ...} containing
                     the a vector of comoving distances ("chi") together with
                     the number counts of source galaxies at these distances
                     or redshifts, respectively.

                     The precision parameter prec is an integer which defines
                     at how many redshifts the matter power spectrum is emula-
                     ted for integration in the Limber equation. The relation
                     between the number of redshifts len_redshifts and the
                     parameter prec is given by

                                     len_redshifts = 2^prec + 1

                     REMARK: The 'chi' vector in the sourcedist dictionary
                     should span the range from 0 to chi(z=5).

        Input type:  emu_pars_dict - dictionary
                     sourcedist - function object
                     prec - int

        Ouput type:  dictionary of the form {'l': ..., 'Cl': ...}
        """
        def limber_prefactor(emu_pars_dict):
            """
            * NESTED FUNCTION * (defined inside get_pconv)

            Signature:          LimberPrefactor(emu_pars_dict)

            Description:        Computes the cosmology dependend prefactor
                                of the integral in the Limber approximation.

            Input type:         emu_pars_dict - dictionary

            Output type:        float
            """
            c_inv = 1./_bg.SPEED_OF_LIGHT_IN_KILOMETERS_PER_SECOND
            h = emu_pars_dict['h']
            H0 = 100*h
            Om_m = emu_pars_dict['om_m']/(h*h)

            # compute prefactor of limber equation integral
            # --> in units of [Mpc^4]
            prefac = 2.25 * Om_m*Om_m * H0*H0*H0*H0 * c_inv*c_inv*c_inv*c_inv

            # convert prefactor units to [(Mpc/h)^4]
            prefac *= h*h*h*h

            return prefac

        def eval_lensingefficiency(emu_pars_dict, sd_func, z_vec):
            """
            * NESTED FUNCTION * (defined inside get_pconv)

            Signature:  eval_lensingefficiency(emu_pars_dict, sd_func, z_vec)

            Descrition: Evaluates the lensing efficiency function for a given
                        cosmology and source galaxy distribution function at
                        a number of different comoving distances (corresponding
                        to different redshifts).

            Input types: emu_pars_dict - dictionary
                         sd_func - function object
                         z_vec - numpy.ndarray

            Output types: numpy.ndarray
                          numpy.ndarray
            """
            chi_lim = _bg.dist_comov(emu_pars_dict, 1e-12, 5.0) # in Mpc/h
            chi_vec = _bg.dist_comov(emu_pars_dict,
                                     _np.zeros_like(z_vec),
                                     z_vec)

            source_dict = {'chi': chi_vec, 'n': sd_func(z_vec)}

            lenseff = _lens.lens_efficiency(source_dict, chi_vec, chi_lim)

            return (chi_vec, lenseff)

        def get_pnonlin_of_l_and_z(emu_pars_dict, z_vec, l_vec):
            """
            * NESTED FUNCTION * (defined inside get_pconv)

            Signature:          get_P_of_l_and_z(emu_pars_dict, z_vec, l_vec)

            Description:        This function computes the different k
                                ranges for the given cosmology at all
                                different redshifts and computes the power
                                spectrum for these redshifts at these k modes.

            Input types:        emu_pars_dict - dictionary
                                z_vec - numpy.ndarray
                                l_vec - numpy.ndarray

            Output type:        numpy.ndarray
            """
            P = get_pnonlin(emu_pars_dict, z_vec) # call EuclidEmulator

            pnonlin_array = []
            for i, z in enumerate(z_vec):
                func = _CubicSpline(_np.log10(P['k']),
                                    _np.log10(P['P_nonlin']['z'+str(i)]))
                k = _cc.l_to_k(emu_pars_dict, l_vec, z) # needs to be called
                                                        # inside loop because
                                                        # different z-values
                                                        # lead to different
                                                        # results

                # evaluate the interpolating function of Pnl for all k in the
                # range allowed by EuclidEmulator (this range is given by
                # P['k'])
                pmatternl = [10.0**func(_np.log10(kk)) for kk in k
                             if kk >= P['k'].min() and kk <= P['k'].max()]
                # for k values below the lower bound or above the upper bound
                # of this k range, set the contributions to 0.0
                p_toosmall = [0.0 for kk in k if kk < P['k'].min()]
                p_toobig = [0.0 for kk in k if kk > P['k'].max()]

                pnonlin_array.append(_np.array(p_toosmall +
                                               pmatternl +
                                               p_toobig))

            # get the non-linear matter power spectrum in units of [(Mpc/h)^3]
            return _np.asarray(pnonlin_array).transpose()

        # =====================================
        if _Class.__module__ not in _sys.modules:
            print("You have not imported neither classee nor classy.\n \
                   Emulating convergence power spectrum is hence not\n \
                   possible.")
            return None

        z_vec = _np.logspace(_np.log10(5e-2),
                             _np.log10(4.999999),
                             2**prec + 1)
        a_inv_vec = 1.+z_vec

        len_l_vec = int(1e4)
        l_vec = _np.linspace(1e1, 2e3, len_l_vec)

        # get the lensing efficiency
        chi_vec, q_vec = eval_lensingefficiency(emu_pars_dict,
                                                sourcedist_func,
                                                z_vec)
        assert len(z_vec) == len(chi_vec)


        # get the matter power spectrum at the correct k modes
        pnonlin_array = get_pnonlin_of_l_and_z(emu_pars_dict, z_vec, l_vec)

        assert all([len(chi_vec) == len(pnonlin_array[i])
                    for i in range(len(pnonlin_array))])

        # Compose integrand
        integrand = _np.array([q_vec*q_vec * a_inv_vec*a_inv_vec *
                               pnonlin_array[i] for i in range(len_l_vec)])

        assert integrand.shape == pnonlin_array.shape

        # perform integral (bear in mind that only the product of the
        # integral and the prefac is truely dimensionless) and return
        return {'l': l_vec,
                'Cl': limber_prefactor(emu_pars_dict)*
                      _romb(integrand, chi_vec[1]-chi_vec[0])}
