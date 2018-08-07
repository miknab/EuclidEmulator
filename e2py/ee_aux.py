"""
ee_aux.py

EuclidEmulator submodule for auxiliary functions.
"""
import sys
import contextlib
import numpy as np

from functools import wraps, partial
from scipy.integrate import quad

# Auxiliary functions for formatting input such that CLASS understands
@contextlib.contextmanager
def disable_numpy_summarization():
    # Author: Jeppe Mosgaard Dakin
    """
    Signature:   disable_numpy_summarization()

    Description: Allows numpy array to have arbitrary length without summarizing
                 its entries.

    Input type:  None

    Output type: None
    """
    threshold = np.get_printoptions()['threshold']
    np.set_printoptions(threshold=np.inf)
    try:
        yield
    finally:
        np.set_printoptions(threshold=threshold)

def stringify_arr(arr):
    # Author: Jeppe Mosgaard Dakin
    """
    Signature:   stringify_arr(arr)

    Description: Takes an array and returns all entries as a single string.

    Input type:  array

    Output type: tuple (array, string)
    """
    with disable_numpy_summarization():
        arr_str = np.array2string(
            arr,
            max_line_width=np.inf,
            formatter={'float': lambda f: '{:.8e}'.format(f)},
            separator=',',
        ).strip('[]')
    return np.fromstring(arr_str, sep=','), arr_str

def print_cosmology(emu_pars_dict):
    # Author: Mischa Knabenhans
    h = emu_pars_dict['h']
    omega_m = emu_pars_dict['om_m']
    omega_rad = 4.183709411969527e-5; # corresponds to 2.755 K Tcmb

    Om_m = omega_m/(h*h)
    Om_rad = omega_rad/(h*h)
    Om_DE = 1.0-Om_m-Om_rad

    print "#"
    print "# Cosmology:"
    print "# dOmega0: ", Om_m 
    print "# dOmegaRad: ", Om_rad 
    print "# dOmegaDE: ", Om_DE

def progress(count, total, status=''):
    # Source: online
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

def normalize(ave_range=[0, np.inf]):
    # Author: Rongchuan Zhao 
    def decorate(func):
        @wraps(func)
        def norm_func(self, *args):
            # reads in some function and returns its 
            # area-normalized version (as a function object)
            inst_func = partial(func, self)
            ToT = quad(inst_func, ave_range[0], ave_range[1], args=args[1:] )[0]
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
