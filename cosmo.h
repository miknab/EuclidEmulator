#ifndef COSMO_HINCLUDED
#define COSMO_HINCLUDED

#define USE_GSL_COSMO
#ifdef USE_GSL_COSMO
#include <gsl/gsl_integration.h>
#endif

    struct csmVariables {
	int bComove;
	double dHubble0;
	double dOmega0;
	double dLambda;
	double dOmegaRad;
	double dOmegab;
	double dOmegaDE;
	double w0;
	double wa;
	double dSigma8;
	double dNormalization;  /* either sigma8 or normalization must be non-zero */
	double dSpectral;
	};


typedef struct csmContext {
    struct csmVariables val;

#ifdef USE_GSL_COSMO
    gsl_integration_workspace *W;
#endif
    } * CSM;
#ifdef __cplusplus
extern "C" {
#endif
    void csmInitialize(CSM *pcsm);
    void csmFinish(CSM csm);
    double csmRadMatEquivalence(CSM csm);

    static inline double csmExp2Hub(CSM csm, double dExp) {
        double dOmegaCurve = 1.0 - csm->val.dOmega0 - csm->val.dLambda - csm->val.dOmegaDE - csm->val.dOmegaRad;

        assert(dExp > 0.0);
        return csm->val.dHubble0
               *sqrt(csm->val.dOmega0*dExp
                     + dOmegaCurve*dExp*dExp
                     + csm->val.dOmegaRad
                     + csm->val.dOmegaDE*pow(dExp,1.0 - 3.0*(csm->val.w0 + csm->val.wa))*exp(-3.0*csm->val.wa*(1.0 - dExp))
                     + csm->val.dLambda*dExp*dExp*dExp*dExp)/(dExp*dExp);
        }

    double csmTime2Hub(CSM csm, double dTime);
    double csmExp2Time(CSM csm, double dExp);
    double csmTime2Exp(CSM csm, double dTime);
    double csmComoveDriftInt(CSM csm, double dIExp);
    double csmComoveKickInt(CSM csm, double dIExp);
    double csmComoveDriftFac(CSM csm, double dTime, double dDelta);
    double csmComoveKickFac(CSM csm, double dTime, double dDelta);
    double csmComoveLookbackTime2Exp(CSM csm, double dComoveTime);  
    void csmComoveGrowth(CSM csm, double a, double *D1LPT, double *D2LPT, double *f1LPT, double *f2LPT);
#ifdef __cplusplus
}
#endif

/*
 ** returns the speed of light in simulation units, given
 ** the simulation length unit in h^-1 Mpc.
 */
static inline double dLightSpeedSim(double dMpcUnit) {
    /*
    ** Find the speed of light in simulation units.
    **
    ** c[Mpc/Gyr] = c[cm/s] * Julian Year[s] / pc[cm] * 1000 
    ** c_sim = c[Mpc/Gyr] * (x Gyrs/ 1 sim time) * ( 1 sim length/Boxsize (Mpc))
    ** x = 1/sqrt(4.498*h*h*2.776e-4)
    */
    /*return(8676.85/dMpcUnit);*/

    /* 
    ** Doug's version:
    **
    ** Cosmological coordinates
    ** G     = 4.30172e-9 Mpc/M. (km/s)^2
    ** rho_c = 3 H^2 / (8 pi G)
    ** c     = 299792.458 km/s
    **
    ** c_sim = c[km/s] * sqrt(Lbox / (G * rho_c * Lbox^3))
    **       = c[km/s] * sqrt(8 pi / (3 H^2 Lbox^2) )
    **       = c[km/s] * sqrt(8 pi / 3) / Lbox / H
    **       = c[km/s] * sqrt(8 pi / 3) / Lbox / h / 100
    ** dMpcUnit given in Mpc/h gives:
    **       = 299792.458 * sqrt(8 pi / 3) / 100 / dMpcUnit
    **       = 8677.2079486362706 / dMpcUnit
    */
    return 8677.2079486362706 / dMpcUnit;
    }

#endif
