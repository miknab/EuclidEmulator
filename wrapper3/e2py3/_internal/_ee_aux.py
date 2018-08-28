"""
_ee_aux.py

EuclidEmulator submodule for auxiliary functions.
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
import contextlib as _ctlib
import numpy as _np

from functools import wraps as _wraps
from functools import partial as _partial
from scipy.integrate import quad as _quad

# Auxiliary functions for formatting input such that CLASS understands
@_ctlib.contextmanager
def disable_numpy_summarization():
    # Author: Jeppe Mosgaard Dakin
    """
    Signature:   disable_numpy_summarization()

    Description: Allows numpy array to have arbitrary length without summarizing
                 its entries.

    Input type:  None

    Output type: None
    """
    threshold = _np.get_printoptions()['threshold']
    _np.set_printoptions(threshold=_np.inf)
    try:
        yield
    finally:
        _np.set_printoptions(threshold=threshold)

def stringify_arr(arr):
    # Author: Jeppe Mosgaard Dakin
    """
    Signature:   stringify_arr(arr)

    Description: Takes an array and returns all entries as a single string.

    Input type:  array

    Output type: tuple (array, string)
    """
    with disable_numpy_summarization():
        arr_str = _np.array2string(
            arr,
            max_line_width=_np.inf,
            formatter={'float': lambda f: '{:.8e}'.format(f)},
            separator=',',
        ).strip('[]')
    return _np.fromstring(arr_str, sep=','), arr_str

def print_cosmology(emu_pars_dict):
    # Author: Mischa Knabenhans
    h = emu_pars_dict['h']
    omega_m = emu_pars_dict['om_m']
    omega_rad = 4.183709411969527e-5; # corresponds to 2.755 K Tcmb

    Om_m = omega_m/(h*h)
    Om_rad = omega_rad/(h*h)
    Om_DE = 1.0-Om_m-Om_rad

    print("#")
    print("# Cosmology:")
    print("# dOmega0: ", Om_m) 
    print("# dOmegaRad: ", Om_rad) 
    print("# dOmegaDE: ", Om_DE)

def progress(count, total, status=''):
    # Source: online
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    _sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    _sys.stdout.flush()

def normalize(ave_range=[0, _np.inf]):
    # Author: Rongchuan Zhao 
    def decorate(func):
        @_wraps(func)
        def norm_func(self, *args):
            # reads in some function and returns its 
            # area-normalized version (as a function object)
            inst_func = _partial(func, self)
            ToT = _quad(inst_func, ave_range[0], ave_range[1], args=args[1:] )[0]
            return func(self, *args)/ToT
        return norm_func
    return decorate

class Function(object):
    # Author: Rongchuan Zhao
    def __init__(self, func):
        self.func = func
    
    def __add__(self, other):
        return Function(lambda x: self.func(x) + other.func(x))
    
    def __sub__(self, other):
        return Function(lambda x: self.func(x) - other.func(x))
    
    def __mul__(self, other):
        return Function(lambda x: self.func(x) * other.func(x))
    
    def __div__(self, other):
        return Function(lambda x: self.func(x) / other.func(x))
    
    def __pow__(self, index):
        return Function(lambda x: self.func(x)**index)
    
    def __call__(self, *arg, **kwarg):
        return self.func(*arg, **kwarg)
