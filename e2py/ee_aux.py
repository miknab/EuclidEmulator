"""
ee_aux.py

EuclidEmulator submodule for auxiliary functions.
"""
import sys
import contextlib
import numpy as np

# Auxiliary functions for formatting input such that CLASS understands
@contextlib.contextmanager
def disable_numpy_summarization():
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
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()
