"""
ee_aux.py

EuclidEmulator submodule for auxiliary functions.
"""

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
            formatter={'float': lambda f: '{:.16e}'.format(f)},
            separator=',',
        ).strip('[]')
    return np.fromstring(arr_str, sep=','), arr_str
