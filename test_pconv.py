import e2py
import numpy as np
import matplotlib.pyplot as plt

MyCosmo = {'om_b': 0.0219961, 'om_m': 0.1431991, 'n_s': 0.96, 'h': 0.67, 'w_0': -1.0, 'sigma_8': 0.83}

nGal = e2py.sgaldist()

pkappa = e2py.get_pconv(MyCosmo, nGal)

l = pkappa['l']
Cl = pkappa['Cl']

y = l*(l+1)*Cl/(2*np.pi)

Fig, ax = plt.subplots()
ax.loglog(l,y)
ax.set_xlabel(r"$\ell$ [1]")
ax.set_ylabel(r"$\frac{\ell(\ell+1)P_{\kappa}(\ell)}{2\pi}$ [1]")
plt.show()
