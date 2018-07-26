import e2py
#import matplotlib.pyplot as plt
import numpy as np

csm = {'om_b': 0.0219961,
       'om_m': 0.1431991,
       'n_s': 0.96,
       'h': 0.67,
       'w_0': -1.0,
       'sigma_8': 0.83}
z=0.0

pnl_dict = e2py.get_pnonlin(csm,z)

kvec = pnl_dict['k']
P_nonlin = pnl_dict['P_nonlin']
P_lin = pnl_dict['P_lin']
Boost = pnl_dict['B']

#np.savetxt("./DataOutput/EmuData.dat", np.c_[kvec, P_nonlin[0], P_lin[0], Boost])

#Fig, axs = plt.subplots(3,1,sharex=True)

#ax = axs[0]
#ax.loglog(kvec,P_lin[0],c="b")
#ax.set_ylabel(r"$P_{\rm lin}(k)\enspace[({\rm Mpc}/h)^3]$")
#
#ax = axs[1]
#ax.axhline(y=1.0,ls=":",c="k")
#ax.loglog(kvec,Boost,c="b")
#ax.set_ylabel(r"$B(k)\enspace[1]$")

#ax = axs[2]
#ax.loglog(kvec,P_nonlin[0],c="b")
#ax.set_ylabel(r"$P_{\rm nl}(k)\enspace[({\rm Mpc}/h)^3]$")
#ax.set_xlabel(r"$k\enspace[h/{\rm Mpc}]$")
#plt.show()
