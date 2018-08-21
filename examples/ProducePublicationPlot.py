import numpy as np
import matplotlib.pyplot as plt
import e2py

from classy import Class

csm = {'om_b': 0.0219961,
       'om_m': 0.1431991,
       'n_s': 0.96,
       'h': 0.67,
       'w_0': -1.0,
       'sigma_8': 0.83}

h = csm['h']
#zvec = np.array([0.0,0.5,1.0,2.0])
Pnl = e2py.get_pnonlin(csm,0.5)
kvec = Pnl['k']
kshape = kvec.shape

ClassyPars = e2py.emu_to_class(csm)
ClassyPars['Omega_Lambda']=0.0
ClassyPars['output']='mPk'
ClassyPars['non linear']='Halofit'
ClassyPars['format']='camb'
ClassyPars['P_k_max_h/Mpc']=10.
ClassyPars['k_per_decade_for_pk']=300.
ClassyPars['z_pk']=0.5#'0.0','0.5','1.0','2.0' 

cosmo=Class()
cosmo.set(ClassyPars)
cosmo.compute()

pHF = np.array([cosmo.pk(k*h,0.5)*h*h*h for k in kvec]).reshape(kshape)

Fig, axs = plt.subplots(2,1, sharex=True)
ax = axs[0]
ax.loglog(kvec,Pnl['P_lin'], c='gray', label = r"$P_\rm{lin}^\rm{CLASS}$")
ax.loglog(kvec,Pnl['P_nonlin'], c='blue', label = r"$P_\rm{nl}^\rm{EE}=P_\rm{lin}^\rm{CLASS} * B$")
ax.grid(True)
ax.set_ylabel(r"P(k,z=0.5) [$(\rm{Mpc}/h)^3$]")
ax.set_xlim([0.01,5])

ax = axs[1]
ax.axhline(y=0, c="black", ls=":")
ax.axhspan(-1,1,color="gray", alpha=0.5)
ax.semilogx(kvec, 100*(Pnl['P_nonlin']/pHF-1), c="black")
ax.grid(True)
ax.set_xlabel(r"$k [h/{\rm Mpc}]$")
ax.set_ylabel(r"$\frac{P_{nl}^{EE}(k,z=0.5)-P^{THM}_{nl}(k,z=0.5)}{P^{THM}_{nl}(k,z=0.5)} [(\rm{Mpc}/h)^3]$")
ax.set_xlim([0.01,5])

plt.show()
plt.show()

  
