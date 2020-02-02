"""
ee_cosmoconv.py

EuclidEmulator submodule for converting cosmological quantities.

REMARK:      The geometry of the Universe is fixed to be flat (i.e.
             Omega_curvature = 1) and the radiation energy density
             is set to Om_rad = 4.183709411969527e-5/(h*h). These
             values were assumed in the construction process of 
             EuclidEmulator and hence must be used whenever one is
             working with it.
"""

import numpy as _np

def z_to_a(z):
    """
    Signature:   z_to_a(z)

    Description: Converts a redshift z into the corresponding
                 scale factor a.

    Input type:  float or numpy.ndarray

    Output type: float or numpy.ndarray
    """
    return 1.0/(z+1.0)

def a_to_z(a):
    """
    Signature:   a_to_z(a)

    Description: Converts a scale factor a into the corresponding
                 redshift z.

    Input type:  float or numpy.ndarray

    Output type: float or numpy.ndarray
    """

    return (1.0/a)-1.0

def sigma8_to_as(sigma8):
    """
    Signature:    sigma8_to_as(sigma8)

    Description:  Converts sigma_8 into A_s based on renormalization.

    Input type:   float

    Output type:  float


    Remarks:      This conversion is not generic. It is the one used
                  in the construction process of Euclid emulator and
                  has to be used whenever a conversion from sigma8
                  to As is performed in the context of working with
                  EuclidEmulator.

    Related:      as_to_sigma8
    """
    # Check if input is given:
    if sigma8 is None:
        raise ValueError("Value of sigma8 is 'None', i.e. sigma8 has \
                          no assigned value.")

    # Check if input is a real number
    sigma8_is_complex = isinstance(sigma8, complex)
    sigma8_is_array = isinstance(sigma8, _np.ndarray)
    complex_entry = False

    if sigma8_is_array:
        complex_entry = any([isinstance(val, _np.complex128) for val in sigma8])

    if sigma8_is_complex or complex_entry:
        raise TypeError("sigma8 must be a real number.")

    a_s = 2.215e-9*(sigma8*sigma8)/(0.8496*0.8496)
    return a_s

def as_to_sigma8(a_s):
    """
    Signature:    as_to_sigma8(a_s)

    Description:  Converts As into sigma8 based on renormalization.

    Input type:   float

    Output type:  float

    Remarks:      This conversion is not generic. It is the one used
                  in the construction process of Euclid emulator and
                  has to be used whenever a conversion from As to
                  sigma8 is performed in the context of working with
                  EuclidEmulator.

    Related:      sigma8_to_as
    """
    sigma_8 = _np.sqrt(0.8496*0.8496*a_s/2.215e-9)
    return sigma_8

def a_to_hubble(emu_pars_dict, a):
    """
    Signature:   a_to_hubble(emu_pars,dict, a, H0, Om_m, Om_DE, Om_curve, w_0, w_a)

    Description: Computes the Hubble parameter at a given scale factor a
                 for a given cosmology.

    Input type:  type(a) = float or numpy.ndarray
                 all other input parameters are of type float

    Output type: type(Hubble) = float or numpy.ndarray

    REMARK:      The geometry of the Universe is fixed to be flat (i.e.
                 Omega_curvature = 1) and the radiation energy density
                 is set to Om_rad = 4.183709411969527e-5/(h*h). These
                 values were assumed in the construction process of 
                 EuclidEmulator and hence must be used whenever one is
                 working with it.
    """
    a_inv = 1.0/a
    
    h = emu_pars_dict['h']
    Om_m = emu_pars_dict['om_m']/(h*h)
    H0 = 100.0 * h # * (km/s)/Mpc
    w_0 = emu_pars_dict['w_0']

    # We explicitly set the radiation energy density. By hardcoding this
    # we make sure that the user cannot choose Om_rad values inconsistent with 
    # that used in the construction process of EuclidEmulator.
    # We adopt the same equation for the computation of the radiation energy 
    # density as in CLASS (cf. CLASS code, input.c file, lines 643 & 714)
    Om_rad = 4.183709411969527e-5/(h*h)

    # We infer Om_DE assuming flat geometry (Om_curve = 0.0). By hardcoding this
    # we make sure that the user cannot choose Om_curve values inconsistent with 
    # that used in the construction process of EuclidEmulator.
    Om_curve = 0.0
    Om_DE = 1.0-Om_curve-Om_m-Om_rad
    
    assert (Om_curve == 0.0 and Om_m + Om_rad + Om_DE + Om_curve == 1.0)

    curvature = Om_curve * a_inv * a_inv
    matter = Om_m * a_inv * a_inv * a_inv
    radiation = Om_rad * a_inv * a_inv * a_inv * a_inv
    darkenergy = Om_DE * a_inv**(3.0+3.0*w_0)

    H = H0 * _np.sqrt(curvature + matter + radiation + darkenergy)

    return H # in the standard units of (km/s)/Mpc

def k_to_l(emu_pars_dict, k, z, prec=12):
    """
    Signature:   k_to_l(emu_pars_dict, k, z, prec=12)

    Description: Converts a numpy.array k into a numpy.array l for a given
                 redshift z. Here, k denotes the wave numbers and l the
                 multipole order. The keyword argument "prec" is a precision
                 parameter for the romberg integration that has to be computed
                 in the course of the comoving distance calculation.

    Input types: type(emu_pars_dict) = dict
                 type(k) = numpy.ndarray
                 type(z) = float
                 type(prec) = int

    Ouput type:  numpy.ndarray

    REMARK:      The geometry of the Universe is fixed to be flat (i.e.
                
                                Omega_curvature = 1 

                 and as a result the comoving angular diameter distance
                 equals the comoving distance) and the radiation energy
                 density is set to 

                        Om_rad = 4.183709411969527e-5/(h*h).

                 These values were assumed in the construction process 
                 of EuclidEmulator and hence must be used whenever one
                 is working with it.
    """
    return k * _bg.dist_comov(emu_pars_dict, 1e-13, z, prec)

def l_to_k(emu_pars_dict, l, z, prec=12):
    """
    Signature:   l_to_k(emu_pars_dict, l, z, prec=12)

    Description: Converts a numpy.array l into a numpy.array k for a given
                 redshift z. Here, k denotes the wave numbers and l the
                 multipole order. The keyword argument "prec" is a precision
                 parameter for the romberg integration that has to be computed
                 in the course of the comoving distance calculation.

    Input types: type(emu_pars_dict) = dict
                 type(l) = numpy.ndarray
                 type(z) = float
                 type(prec) = int

    Ouput type:  numpy.ndarray

    REMARK:      The geometry of the Universe is fixed to be flat (i.e.
                
                                Omega_curvature = 1 

                 and as a result the comoving angular diameter distance
                 equals the comoving distance) and the radiation energy
                 density is set to 

                        Om_rad = 4.183709411969527e-5/(h*h).

                 These values were assumed in the construction process 
                 of EuclidEmulator and hence must be used whenever one
                 is working with it.
    """
    return l/_bg.dist_comov(emu_pars_dict, 1e-13, z, prec)
