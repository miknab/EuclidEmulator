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

import cython
cimport cython
from libc.stdlib cimport free
import sys
import numpy as np

#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

C2np = {
    # Booleans
    'bint': np.bool,
    # Integers
    'char'         : np.byte,
    'short'        : np.short,
    'int'          : np.intc,
    'long int'     : np.long,
    'long long int': np.longlong,
    'ptrdiff_t'    : np.intp,
    'Py_ssize_t'   : np.intp,
    # Unsgined integers
    'unsigned char'         : np.ubyte,
    'unsigned short'        : np.ushort,
    'unsigned int'          : np.uintc,
    'unsigned long int'     : np.uint,
    'unsigned long long int': np.ulonglong,
    'size_t'                : np.uintp,
    # Floating-point numbers
    'float'     : np.single,
    'double'    : np.double,
    'long float': np.longfloat,
}

# Include C-library functions
cdef extern from "sys/mman.h":
    munmap (void *__addr, size_t __len); 

# Import the C function
cdef extern from "EuclidEmulator.h": 
    struct FID:
        double *handle
        int size
    FID EucEmu(double *CosmoParams, double *Redshifts, int len_z, double **kVals, int *nkVals, double **Boost, int *nBoost)

# Simple class storing the data passed back from the
# wrapped C function. It also takes care of freeing the data.
cdef class Py_EE_class:
    cdef:
        # Internal data
        double[::1] mv_k
        double[::1] mv_b
        double *ptr0
        double *ptr1
        
        # Variables accessible from Python space
        public object k
        public object boost
    
    def __cinit__(self, double[::1] mv_k, double[::1] mv_b):
        # Store arrays as both memoyviews and pointers
        self.mv_k = mv_k
        self.mv_b = mv_b
        self.ptr0 = &self.mv_k[0]
        self.ptr1 = &self.mv_b[0]

        # Wrap the Cython memoryviews in NumPy arrays
        self.k = np.asarray(self.mv_k)
        self.boost = np.asarray(self.mv_b) 


    # *** NOTE ***  
    # The following __dealloc__ function is not needed as the memory mapped data in the C function
    # is also unmapped inside the C function and all other arrays and memory views are garbage collected
    #

    #def __dealloc__(self):
    #    if self.ptr0 is not NULL:
    #        print("right before free:", self.fsize)    
    #        retcode = munmap(self.ptr0, self.fsize) 
    #        print("retcode: %d" %retcode)
    #    if self.ptr1 is not NULL:
    #        free(self.ptr1)

# Function with both C and Python capablities,
# wrapping the C function.
def emu_boost(cosmo_par_array,redshift_array):
    # Cast input variables to Cython types
    cdef double[::1] redshifts = np.asarray(redshift_array, dtype=C2np['double']) 
    cdef int nz = len(redshifts)
    
    cdef double[::1] cosmo_params = np.asarray(cosmo_par_array, dtype=C2np['double'])
 
    # Define output data containers for C-function
    cdef int nkVals = 0
    cdef int nBoost = 0
    cdef double **Boost
    cdef double **kVals 

    cdef double somedouble = 1.0 
    cdef double* somedouble_ptr = &somedouble
    Boost = &somedouble_ptr
    cdef double hmm = 1.0
    cdef double* hmm_ptr = &hmm
    kVals = &hmm_ptr

    # Call the C function
    EucEmu(&cosmo_params[0],&redshifts[0],nz,kVals, &nkVals, Boost, &nBoost)

    # Wrap the data in a simple class. We have to pass the data
    # as Cython memoryviews, as structs (and pointers) are too
    # low level for them to be allowed as arguments to __cinit__.

    cdef double[::1] mv_k, mv_b, mv_f
    mv_k = <double[:nkVals]>(kVals[0])
    mv_b = <double[:nz*nBoost]>(Boost[0]) 
    
    return Py_EE_class(mv_k, mv_b)
