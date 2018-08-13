import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

"""
PROGRAM:     GetPnonlin.py

SYNOPSIS:    python GetPnonlin.py

DESCRIPTION: This python script reads in three data files containing 
             the linear power spectrum (produced by CLASS), the non-
             linear power spectrum (produced by CLASS and Halofit), 
             and the boost factor (produced by EuclidEmulator). From
             the linear power spectrum and the boost, it computes the
             fully non-linear power spectrum and compares it to the 
             Halofit solution.

             The resulting non-linear power spectrum is stored in an
             output file.

             A plot is generated that allows to visually compare the 
             EuclidEmulator and the Halofit solutions. For the original
             parameter settings (Euclid reference cosmology and redshift
             z=0.5), the resulting plot is part of Fig. 9 in the paper
             by Knabenhans et al. (2018).
      
AUTHOR:      M. Knabenhans

DATE:        May, 2018
"""

# Read in data files
# ===================
CLASS_File = "EuclidRef.pk.dat"                 # contains two columns: k P(k)_linear
THM_File   = "EuclidRef.pk_nl.dat"              # contains two columns: k P(k)_non-linear
EE_File    = "EXAMPLE_EucRefBoost.z0.5.dat"     # contains two columns: k B(k)

k_Class, P_lin = np.loadtxt(CLASS_File, unpack=True)
k_THM, P_THM = np.loadtxt(THM_File, unpack=True)     # THM = Takahashi halo model, cf. Takahashi et al. (2012)
k_EE, B = np.loadtxt(EE_File,unpack=True)            # EE = EuclidEmulator, cf. Knabenhans et al. (2018)

# Interpolation of linear power spectrum
# ======================================
# We use a cubic spline interpolation in log space

func_P_lin   = CubicSpline(np.log10(k_Class),np.log10(P_lin))
interp_P_lin = 10**func_P_lin(np.log10(k_EE))

func_P_THM   = CubicSpline(np.log10(k_THM),np.log10(P_THM))
interp_P_THM = 10**func_P_THM(np.log10(k_EE))

# Compute fully non-linear power spectrum
# =======================================
P_nonlin = interp_P_lin * B

# Compare EuclidEmulator to Takahashi
# =======================================
ratio = 100*(P_nonlin/interp_P_THM-1)

# Store data in output file
# =========================
Output_File = "EXAMPLE_P_nonlinear.z0.5.dat"
np.savetxt(Output_File,np.c_[k_EE,P_nonlin])

# Plot results and compare to Takahashi's halo model extension
# ============================================================
Fig, axs = plt.subplots(2,1,sharex=True)

ax = axs[0]
ax.loglog(k_EE,P_nonlin, color = 'blue', label=r"$P_{\rm nl}^{\rm EE} = P_{\rm lin}^{\rm CLASS} * B$")
ax.loglog(k_EE,interp_P_lin, color='gray', label=r"$P_{\rm lin}^{\rm CLASS}$")
ax.grid(True)
ax.legend(loc="lower left")
ax.set_xlim(0.01,5)
ax.set_ylabel(r"P(k,z=0.5) [$({\rm Mpc}/h)^{3}$]")

ax = axs[1]
ax.axhline(y=0.0, ls=":", color="black")
ax.axhspan(-1,1, color="gray", alpha=0.5)
ax.semilogx(k_EE,ratio, color="black")
ax.grid(True)
ax.set_xlabel(r"k [$h/{\rm Mpc}$]")
ax.set_ylabel(r"$\frac{P_{\rm nl}^{\rm EE}(k,z=0.5)-P_{\rm nl}^{\rm THM}(k,z=0.5)}{P_{\rm nl}^{\rm THM}(k,z=0.5)}$ [1]")
plt.tight_layout()
plt.show()
