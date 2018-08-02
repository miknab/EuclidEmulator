import e2py
import numpy as np

chivec=np.linspace(0,44.5,1000)
nvec = chivec**2.0
sdist = {'chi': chivec, 'n': nvec}
MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}
e2py.get_pconv(MyCosmo, sdist, prec=5)
