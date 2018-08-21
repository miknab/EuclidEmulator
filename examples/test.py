import e2py
import matplotlib.pyplot as plt
import numpy as np
import pylab as plb
from scipy.interpolate import CubicSpline

# Specify cosmology and redshifts at which the non-linear
# power spectrum shall be emulated
csm = {'om_b': 0.0219961,
       'om_m': 0.1431991,
       'n_s': 0.96,
       'h': 0.67,
       'w_0': -1.0,
       'sigma_8': 0.83}

z=np.array([0.0,1.0,5.0])

# Emulation step: the one next line does all the magic
pnl_dict = e2py.get_pnonlin(csm,z)

# We rename some variables to avoid long variables names
# later on (this step is clearly optional)
kvec = pnl_dict['k']
P_nonlin = pnl_dict['P_nonlin']
P_lin = pnl_dict['P_lin']
Boost = pnl_dict['B']

# Save the data
for idx in P_nonlin.keys():
    np.savetxt("./DataOutput/EmuData."+idx+".dat", np.c_[kvec, P_nonlin[idx], P_lin[idx], Boost[idx]])

# Simulation data
ksim, Psim, __ = np.loadtxt("AVG_EuclidReference.wRad.00100.pk", unpack=True)
fsim = CubicSpline(np.log10(ksim), np.log10(Psim))
Psim_intp = 10.0**fsim(np.log10(kvec))

# Plot the data
Fig, axs = plt.subplots(4,1,sharex=True)

ax = axs[0]
ax.loglog(kvec,P_lin['z0'],c="b",label='z='+str(z[0]))
ax.loglog(kvec,P_lin['z1'],c="r",label='z='+str(z[1]))
ax.loglog(kvec,P_lin['z2'],c="g",label='z='+str(z[2]))
ax.legend(loc='lower left')
ax.set_ylabel(r"$P_{\rm lin}(k)\enspace[({\rm Mpc}/h)^3]$")

ax = axs[1]
ax.axhline(y=1.0,ls=":",c="k")
ax.loglog(kvec,Boost['z0'],c="b",label='z='+str(z[0]))
ax.loglog(kvec,Boost['z1'],c="r",label='z='+str(z[1]))
ax.loglog(kvec,Boost['z2'],c="g",label='z='+str(z[2]))
ax.legend(loc='upper left')
ax.set_ylabel(r"$B(k)\enspace[1]$")

ax = axs[2]
ax.loglog(kvec,P_nonlin['z0'],c="b",label='z='+str(z[0]))
ax.loglog(kvec,P_nonlin['z1'],c="r",label='z='+str(z[1]))
ax.loglog(kvec,P_nonlin['z2'],c="g",label='z='+str(z[2]))
ax.loglog(kvec,Psim_intp,c='k', ls="--", label="z=0.0 (sim/ref)")
ax.legend(loc='lower left')
ax.set_ylabel(r"$P_{\rm nl}(k)\enspace[({\rm Mpc}/h)^3]$")
#ax.set_xlabel(r"$k\enspace[h/{\rm Mpc}]$")

ax = axs[3]
ax.axhspan(-1, 1, alpha=0.5, color='grey')
ax.axhline(y=0.0, ls=":", c='k')
ax.semilogx(kvec, 100*(P_nonlin['z0']/Psim_intp-1), c='b', ls="-")
ax.set_ylabel(r"$\varepsilon_{rel}\;[\%]$")
ax.set_xlabel(r"$k\enspace[h/{\rm Mpc}]$")
ax.set_xlim([0.005,5.35])

plt.show()


