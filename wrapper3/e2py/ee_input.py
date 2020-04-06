"""
ee_input.py

EuclidEmulator submodule containing functions related to argument parsing.
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
import pandas as _pd
from e2py._internal import _ee_cosmoconv as _cc
######################################################
#################### Check input #####################
######################################################

def check_param_range(par_dict, csm_index=0):
    """
    Checks if all parameters in the cosmology dictionary 'par_dict'
    (with index 'csm_index') passed to this function obey the limits
    set by the emulator.
    """
    om_b_range = [0.0217, 0.0233]
    om_m_range = [0.1326, 0.1526]
    n_s_range = [0.9345, 0.9965]
    h_range = [0.6251, 0.7211]
    w_0_range = [-1.250, -0.750]
    sigma_8_range = [0.7684, 0.8614]

    om_b_not_in_range = om_b_range[0] > par_dict['om_b'] or\
                        om_b_range[1] < par_dict['om_b']

    om_m_not_in_range = om_m_range[0] > par_dict['om_m'] or\
                        om_m_range[1] < par_dict['om_m']

    n_s_not_in_range = n_s_range[0] > par_dict['n_s'] or\
                       n_s_range[1] < par_dict['n_s']

    h_not_in_range = h_range[0] > par_dict['h'] or\
                     h_range[1] < par_dict['h']

    w_0_not_in_range = w_0_range[0] > par_dict['w_0'] or\
                       w_0_range[1] < par_dict['w_0']

    sigma_8_not_in_range = sigma_8_range[0] > par_dict['sigma_8'] or\
                           sigma_8_range[1] < par_dict['sigma_8']

    if om_b_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          om_b is set to %f, but should be in the interval \
                          [0.0217, 0.0233]."
                         %(csm_index, par_dict['om_b']))

    if om_m_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          om_m is set to %f, but should be in the interval \
                          [0.1326, 0.1526]."
                         %(csm_index, par_dict['om_m']))

    if n_s_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          n_s is set to %f, but should be in the interval \
                          [0.9345, 0.9965]."
                         %(csm_index, par_dict['n_s']))

    if h_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          h is set to %f, but should be in the interval \
                          [0.6251, 0.9965]."
                         %(csm_index, par_dict['h']))

    if w_0_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          w_0 is set to %f, but should be in the interval \
                          [-1.250, -0.750]."
                         %(csm_index, par_dict['w_0']))

    if sigma_8_not_in_range:
        raise ValueError("Parameter range violation in cosmology %d: \
                          sigma_8 is set to %f, but should be in the interval \
                          [0.7684, 0.8614]."
                         %(csm_index, par_dict['sigma_8']))

######################################################
#################### Argument parsing ################
######################################################

def read_parfile(filename, sep=","):
    """
    Signature:   read_parfile(filename, sep=",")

    Description: Reads in a parameter file and returns a dictionary or a
                 list of dictionaires (if multiple cosmologies are specified)
                 with the parameters.

    Input type:  string (path to or name of parameter file)

    Output type: list of python dictionaries
    """
    list_of_cosmologies = []
    parameters = _pd.read_csv(filename,
                              delimiter='\s*'+sep+'\s*',
                              engine='python')

    for indx, row in parameters.iterrows():
        # Convert pandas.Series to dictionary (= required output type)
        cosmo_dict = row.to_dict()
        # Check if parameter ranges are obeyed
        check_param_range(cosmo_dict, indx)
        # Append to list
        list_of_cosmologies.append(cosmo_dict)

    return list_of_cosmologies

########################################################
################### Format conversion ##################
########################################################

def emu_to_class(emu_pars_dict):
    """
    Signature:    emu_to_class(emu_pars_dict)

    Description:  Converts the set of parameters accepted by
                  EuclidEmulator into a set of parameters accepted
                  by CLASS and CAMB.

    Input type:   python dictionary

    Output type:  python dictionary

    Remark:       The expected keys are:
                  'om_b'    (lowercase baryon density parameter)
                  'om_m'    (lowercase total matter density parameter)
                  'n_s'     (spectral index)
                  'h'       (Hubble parameter)
                  'w_0'     (dark energy equation of state parameter)
                  'sigma_8' (overdensity fluctuation variance).

    Related:      class_to_emu
    """
    if not isinstance(emu_pars_dict, dict):
        print("The cosmological parameters must be passed as a \
               python dictionary.\n")
        _sys.exit()

    om_b = emu_pars_dict['om_b']
    om_m = emu_pars_dict['om_m']
    n_s = emu_pars_dict['n_s']
    h = emu_pars_dict['h']
    w_0 = emu_pars_dict['w_0']
    sigma_8 = emu_pars_dict['sigma_8']

    om_cdm = om_m - om_b

    class_pars_dict = {'omega_b': om_b,
                       'omega_cdm': om_cdm,
                       'n_s': n_s,
                       'h': h,
                       'w0_fld': w_0,
                       'sigma8': sigma_8
                      }

    return class_pars_dict

def class_to_emu(class_pars_dict):
    """
    Signature:    class_to_emu(class_pars_dict)

    Description:  Converts the set of parameters accepted by CLASS
                  and CAMB into a set of parameters accepted by
                  EuclidEmulator.

    Input type:   python dictionary

    Output type:  python dictionary

    Remark:       The expected keys are:
                  'omega_b'   (lowercase baryon density parameter)
                  'omega_cdm' (lowercase cold dark matter density parameter)
                  'n_s'       (spectral index)
                  'h'         (Hubble parameter)
                  'w0_fld'    (dark energy equation of state parameter)
                  'sigma8     (amplitude of power spectrum/variance of density field).

    Related:      class_to_emu
    """
    if not isinstance(class_pars_dict, dict):
        print("The cosmological parameters must be passed as a \
               python dictionary.\n")
        _sys.exit()

    om_b = class_pars_dict['omega_b']
    om_cdm = class_pars_dict['omega_cdm']
    n_s = class_pars_dict['n_s']
    h = class_pars_dict['h']
    w_0 = class_pars_dict['w0_fld']
    sigma_8 = class_pars_dict['sigma8']

    om_m = om_b + om_cdm

    emu_pars_dict = {'om_b': om_b,
                     'om_m': om_m,
                     'n_s': n_s,
                     'h': h,
                     'w_0': w_0,
                     'sigma_8': sigma_8}

    return emu_pars_dict
