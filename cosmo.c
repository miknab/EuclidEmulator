#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif
#include "cosmo.h"

//This version of the code requires GSL
#ifndef USE_GSL_COSMO
#error USE_GSL_COSMO must be defined!
#else
#define LIMIT 1000
#endif

/*
 * Cosmological module for PKDGRAV.
 * N.B.  This code is being shared with skid and the I.C. generator.
 */

void csmInitialize(CSM *pcsm) {
    CSM csm;

    csm = (CSM) malloc(sizeof(struct csmContext));
    assert(csm != NULL);

    csm->val.dHubble0 = 0.0;
    csm->val.dOmega0 = 0.0;
    csm->val.dLambda = 0.0;
    csm->val.dOmegaDE = 0.0;
    csm->val.w0 = 0.0;
    csm->val.wa = 0.0;
    csm->val.dOmegaRad = 0.0;
    csm->val.dOmegab = 0.0;
    csm->val.bComove = 0;
    csm->W = gsl_integration_workspace_alloc(LIMIT);
    *pcsm = csm;
    }

void csmFinish(CSM csm) {
    gsl_integration_workspace_free(csm->W);
    free(csm);
    }

#define EPSCOSMO 1e-7

/*
 * ** by MK: Computes the scale factor a at radiation-matter equivalence.
 * */
double csmRadMatEquivalence(CSM csm){
    return csm->val.dOmegaRad/csm->val.dOmega0;
}

double csmTime2Hub(CSM csm,double dTime) {
    double a = csmTime2Exp(csm,dTime);

    assert(a > 0.0);
    return csmExp2Hub(csm, a);
    }

static double Exp2Time_integrand(double ak, void * params) {
    CSM csm = (CSM)params;

    double dExp = pow(ak,2.0/3.0);
    assert (dExp > 0.0);

    return 2.0/(3.0*ak*csmExp2Hub(csm,dExp));
    }

static double Exp2TimeIntegrate(CSM csm,double dExp) {
    gsl_function F;
    F.function = &Exp2Time_integrand;
    F.params = csm;
    double result,error;
    gsl_integration_qag(&F, 0.0, pow(dExp, 1.5),
	0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
    //printf("a=%g,\t result of Exp2TimeIntegrate = %g\n", dExp, result);
    return result;
    }

double csmExp2Time(CSM csm,double dExp) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,eta;

    if (!csm->val.bComove) {
	/*
	 ** Invalid call!
	 */
	assert(0);
	}
    
    if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
	if (dOmega0 == 1.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(2.0/(3.0*dHubble0)*pow(dExp,1.5));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		B = 1.0/sqrt(dOmega0);
		eta = acos(1.0-dExp);
		return(B*(eta-sin(eta)));
		}
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
	    A = 0.5*dOmega0/(dOmega0-1.0);
	    B = A*a0;
	    eta = acos(1.0-dExp/A);
	    return(B*(eta-sin(eta)));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta = acosh(dExp/A+1.0);
	    return(B*(sinh(eta)-eta));
	    }
	else if (dOmega0 == 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(dExp/dHubble0);
	    }
	else {
	    /*
	     ** Bad value.
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
	return Exp2TimeIntegrate(csm,dExp);
	}
    }

#define MAX_ITER 100

double csmTime2Exp(CSM csm,double dTime) {
    double al=0,ah=1,a0,a1=1,at,a;
    double th,f,f1,h,ho;
    int j;

    if (!csm->val.bComove) return(1.0);
    else {
	assert(dTime > 0);
	th = csmExp2Time(csm,ah);
	/*
	** Search for upper bracket if needed.
	*/
	while (dTime > th) {
	    a0 = a1;
	    a1 = ah;
	    ah = a1+a0;
	    th = csmExp2Time(csm,ah);
	    }
	a = 0.5*(al+ah);
	ho = ah-al;
	h = ho;

	f = dTime - Exp2TimeIntegrate(csm,a);
	f1 = 1/(a*csmExp2Hub(csm,a));
	for (j=0;j<MAX_ITER;++j) {
	    if (a+f/f1 < al || a+f/f1 > ah || fabs(2*f) > fabs(ho*f1)) {
		/*
		** Bisection Step.
		*/
		ho = h;
		h = 0.5*(ah-al);
		a = al+h;
		/*
				printf("bisect al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == al) return a;
		}
	    else {
		/*
		** Newton Step.
		*/
		ho = h;
		h = f/f1;
		at = a;
		a += h;
		/*
				printf("newton al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == at) return a;
		}
	    if (fabs(h) < EPSCOSMO) {
		/*
				printf("converged al:%.14g ah:%.14g a:%.14g t:%.14g == %.14g\n",
				       al,ah,a,dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,0.0,pow(a,1.5),EPSCOSMO*1e-1),
				       dTime);
		*/
		return a;
		}

	    f = dTime - Exp2TimeIntegrate(csm,a);
	    f1 = 1/(a*csmExp2Hub(csm,a));
	    if (f < 0) ah = a;
	    else al = a;
	    }
	assert(0);
	}
    return 0.0; /* We never reach here, but this keeps the compiler happy */
    }


double csmComoveDriftInt(CSM csm, double dIExp) {
    return -dIExp/(csmExp2Hub(csm, 1.0/dIExp));
    }
static double ComoveDrift_integrand(double diExp, void * params) {
    return csmComoveDriftInt(params,diExp);
    }

/*
 ** Make the substitution y = 1/a to integrate da/(a^2*H(a))
 */
double csmComoveKickInt(CSM csm, double dIExp) {
    return -1.0/(csmExp2Hub(csm, 1.0/dIExp));
    }

static double ComoveKick_integrand(double diExp, void * params) {
    return csmComoveKickInt(params,diExp);
    }

/*
 ** This function integrates the time dependence of the "drift"-Hamiltonian.
 */
double csmComoveDriftFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
 
    else if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
	a1 = csmTime2Exp(csm,dTime);
	a2 = csmTime2Exp(csm,dTime+dDelta);
	if (dOmega0 == 1.0) {
	    return((2.0/dHubble0)*(1.0/sqrt(a1) - 1.0/sqrt(a2)));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		A = 1.0;
		B = 1.0/sqrt(dOmega0);
		}
	    else {
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		}
	    eta1 = acos(1.0-a1/A);
	    eta2 = acos(1.0-a2/A);
	    return(B/A/A*(1.0/tan(0.5*eta1) - 1.0/tan(0.5*eta2)));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta1 = acosh(a1/A+1.0);
	    eta2 = acosh(a2/A+1.0);
	    return(B/A/A*(1.0/tanh(0.5*eta1) - 1.0/tanh(0.5*eta2)));
	    }
	else if (dOmega0 == 0.0) {
	    /*
	     ** YOU figure this one out!
	     */
	    assert(0);
	    return(0.0);
	    }
	else {
	    /*
	     ** Bad value?
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
	gsl_function F;
	F.function = &ComoveDrift_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
	}
    }


/*
 ** This function integrates the time dependence of the "kick"-Hamiltonian.
 */
double csmComoveKickFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
    else if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
	a1 = csmTime2Exp(csm,dTime);
	a2 = csmTime2Exp(csm,dTime+dDelta);
	if (dOmega0 == 1.0) {
	    return((2.0/dHubble0)*(sqrt(a2) - sqrt(a1)));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		A = 1.0;
		B = 1.0/sqrt(dOmega0);
		}
	    else {
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		}
	    eta1 = acos(1.0-a1/A);
	    eta2 = acos(1.0-a2/A);
	    return(B/A*(eta2 - eta1));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta1 = acosh(a1/A+1.0);
	    eta2 = acosh(a2/A+1.0);
	    return(B/A*(eta2 - eta1));
	    }
	else if (dOmega0 == 0.0) {
	    /*
	     ** YOU figure this one out!
	     */
	    assert(0);
	    return(0.0);
	    }
	else {
	    /*
	     ** Bad value?
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
	gsl_function F;
	F.function = &ComoveKick_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
	}
    }

double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime) {
    if (!csm->val.bComove) return(1.0);
    else {
	double dExpOld = 0.0;
	double dT0 = csmExp2Time(csm, 1.0);
	double dTime = dT0 - dComoveTime;
	double dExpNew;
	int it = 0;

	if (dTime < EPSCOSMO) dTime = EPSCOSMO;
	dExpNew = csmTime2Exp(csm, dTime);
	/*
	 * Root find with Newton's method.
	 */
	do {
	    double dTimeNew = csmExp2Time(csm, dExpNew);
	    double f = dComoveTime
		       - csmComoveKickFac(csm, dTimeNew, dT0 - dTimeNew);
	    double fprime = -1.0/(dExpNew*dExpNew*csmExp2Hub(csm, dExpNew));
	    dExpOld = dExpNew;
	    dExpNew += f/fprime;
	    it++;
	    assert(it < 20);
	    }
	while (fabs(dExpNew - dExpOld)/dExpNew > EPSCOSMO);
	return dExpNew;
	}
    }

static double RK4_f1(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}

static double RK4_g1(CSM csm, double lna, double D, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * D/csmExp2Hub(csm, a);
}

// This function is in principle redundant as it is exactly the same as RK4_f1
static double RK4_f2(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}

static double RK4_g2(CSM csm, double lna, double D1, double D2, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * (D2 - D1*D1)/csmExp2Hub(csm, a);
}


void csmComoveGrowth(CSM csm, double a, double *D1LPT, double *D2LPT, double *f1LPT, double *f2LPT){
    /*
    ** Variable declarations & initializations
    */
    int nSteps=1000;
    double a_init, lna_init = log(1e-12); // ln(a)=-12 ==> a = e^(-12) ~ 0 
    double stepwidth = (log(a)- lna_init)/nSteps;

    // NOTICE: Storing the following quantities into data structures is by far not optimal (we actually never need the old values after the update).
    double ln_timesteps[nSteps+1];
    // -- 1LPT 
    double D1[nSteps+1]; // 1LPT Growth factor D1(a)
    double G1[nSteps+1]; // G1(a) = dD1(a)/dln(a) *H  ==>  Growth rate: f1(a) = G1/(H*D1) 
    // -- 2LPT
    double D2[nSteps+1]; // 2LPT Growth factor D2(a)
    double G2[nSteps+1]; // G2(a) = dD2(a)/dln(a) *H  ==>  Growth rate: f2(a) = G1/(H*D1) 

    /* 
    ** Set boundary conditions
    */
    a_init = exp(lna_init);
    D1[0] = a_init + 2.0/3.0*csmRadMatEquivalence(csm); 
    G1[0] = csmExp2Hub(csm, a_init)*a_init;

    // This is the analytical approximation
    //double AnApprox = -3. * D1[0]*D1[0]/(7. * pow(csm->val.dOmega0,1./143.));

    D2[0] = -2.0/3.0*csmRadMatEquivalence(csm)*a_init;
    G2[0] = -2.0/3.0*csmRadMatEquivalence(csm)*a_init*csmExp2Hub(csm, a_init);


    //Classical RK4 Solver
    double k0, k1, k2, k3;
    double l0, l1, l2, l3;
    double m0, m1, m2, m3;
    double n0, n1, n2, n3;

    //FILE *fp;
    //fp = fopen("GrowthFactorTable.NewBC.dat","a");

    int i; // running loop variable
    for(i=0;i<=nSteps;i++){
        ln_timesteps[i] = lna_init + i*stepwidth;
        //fprintf(file, "%.15f, %.5f,%.20f\n", exp(ln_timesteps[i]),1.0/exp(ln_timesteps[i])-1.0, D[i]+ 0.0001977011);
 
        //RK4 step 1
        k0 = stepwidth * RK4_f1(csm, ln_timesteps[i], G1[i]);
        l0 = stepwidth * RK4_g1(csm, ln_timesteps[i], D1[i], G1[i]);

	m0 = stepwidth * RK4_f2(csm, ln_timesteps[i], G2[i]);
        n0 = stepwidth * RK4_g2(csm, ln_timesteps[i], D1[i], D2[i], G2[i]);

	//RK4 step 2
	k1 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth/2.0, G1[i] + l0/2.0); 
	l1 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k0/2.0, G1[i] + l0/2.0);
    
	m1 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth/2.0, G2[i] + n0/2.0);
        n1 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k0/2.0, D2[i] + m0/2.0, G2[i] + n0/2.0);

      	//RK4 step 3
	k2 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth/2.0, G1[i] + l1/2.0);
  	l2 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k1/2.0, G1[i] + l1/2.0);

	m2 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth/2.0, G2[i] + n1/2.0);
        n2 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k1/2.0, D2[i] + m1/2.0, G2[i] + n1/2.0);

	//RK4 step 4
	k3 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth, G1[i] + l2);
	l3 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth, D1[i] + k2, G1[i] + l2);

	m3 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth, G2[i] + n2);
        n3 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth, D1[i] + k2, D2[i] + m2, G2[i] + n2);

	//Update
	D1[i+1] = D1[i] + (k0 + 2*k1 + 2*k2 + k3)/6.0;
	G1[i+1] = G1[i] + (l0 + 2*l1 + 2*l2 + l3)/6.0; 

        D2[i+1] = D2[i] + (m0 + 2*m1 + 2*m2 + m3)/6.0;
        G2[i+1] = G2[i] + (n0 + 2*n1 + 2*n2 + n3)/6.0; 

	//fprintf(fp, "%.20g, %.20g, %.20g\n", exp(lna_init + i*stepwidth), D1[i], D2[i]);
    }
       
    //fclose(fp);

    *D1LPT = D1[nSteps];
    *f1LPT = G1[nSteps]/(csmExp2Hub(csm,a) * *D1LPT);

    *D2LPT = D2[nSteps];
    *f2LPT = G2[nSteps]/(csmExp2Hub(csm,a) * *D2LPT);
    return;
}

