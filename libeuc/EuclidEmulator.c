#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __APPLE__
    #include <sys/uio.h>
#else
    #include <sys/io.h>
#endif

#include <sys/mman.h>
#include <assert.h>
#include <gsl/gsl_sf_legendre.h>
#include "cosmo.h"
#include "EuclidEmulator.h"

void EucEmu(double *CosmoParams, double *Redshifts, int len_z, double **kVals, int *nkVals, double **Boost, int *nBoost){

    const int nSteps = 100;
    const int nk = 1099;
    const int nz = nSteps+1;
    const double omega_rad = 4.183709411969527e-5; /* corresponds to 2.755 K Tcmb */
    const int bVerbose = 0;
    const int nCoef[11] = {175,82,82,82,88,139,82,82,139,82,40};//{76,76,76,82,82,199,82,82};
    const double min[6] = {0.0215,0.1306,0.9283,0.6155,-1.3,0.7591};
    const double max[6] = {0.0235,0.1546,1.0027,0.7307,-0.7,0.8707};

    CSM csm;

    double *coef[11];
    double *index[11];
    double lamb[11];
    double *pca;
    double *mean[nz];
    double *Pl[6]; /* the legendre polynomials */
    double *kvec;
    double *tmp;
    double *lnboost;
    double *boost;
    double *f;
    off_t size;
    struct stat s;

    printf("data file path = %s\n", PATH_TO_EUCLIDEMULATOR_DATA_FILE);

    int fd = open(PATH_TO_EUCLIDEMULATOR_DATA_FILE"/ee.dat", O_RDONLY);
    int zcounter,i,di,ip,iz,ic,j,l,lmax;
    double z,x,prod,t0,t200,dDelta,t,dti;

    assert(fd > 0); 
    /* Get the size of the file. */
    //printf("Get file size...\n");
    int status = fstat(fd, & s);
    size = s.st_size;

    f = (double *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);

    struct FID retFid = {f,size};

    /*
    ** Impose a little structure on the flat file of doubles.
    ** We are just setting pointers here!
    */
    //printf("Imposing structure...\n");
    di = 0;                             /* counts number of floats in data file */
    
    /* Reading in PCE coefficients */
    for (i=0;i<11;++i) {
      coef[i] = &f[di];
      di += nCoef[i];
    }

    /* Reading in PCE multi-indices */
    for (i=0;i<11;++i) {
	index[i] = &f[di];
	di += 6*nCoef[i];
    }

    /* Reading in principal components */
    pca = &f[di];
    di += nk*nz*11;  // skip over the whole "block"

    /* Reading in PCA means */
    for (iz=0;iz<nz;++iz) {
      mean[iz] = &f[di];
      di += nk;
    }

    /* Reading in k vector */
    kvec = malloc(nk*sizeof(double));
    tmp = &f[di];
    for (i = 0; i < nk; i++) {
      kvec[i] = tmp[i];
    }
    di += nk;  

    assert(di == size/sizeof(double));
 
    /*
    ** Determine lmax (= maximal order of Legendre polynomials).
    */
    //printf("Determining lmax...\n");
    lmax = 0;
    for (i=0;i<11;++i) {
      for (ip=0;ip<6;++ip) {
	for (ic=0;ic<nCoef[i];++ic) {
          l = index[i][ic*6 + ip];
	  if (l > lmax) lmax = l;
	}
      }
    }
    /*
    ** Print out the interaction matrix. [optional]
    */
    if (bVerbose) {
      //printf("Printing interaction matrix...\n");
      fprintf(stderr,"lmax = %d\n",lmax);
      for (i=0;i<11;++i) {
	fprintf(stderr,"%1d ",i);
	for (j=0;j<2*nCoef[i]-1;++j) fprintf(stderr,"-");
	fprintf(stderr,"\n");
	for (ic=0;ic<nCoef[i];++ic) {
	    fprintf(stderr,"%.3g ",coef[i][ic]);
	    }
	fprintf(stderr,"\n");
	for (j=0;j<2*nCoef[i]-1;++j) fprintf(stderr,"-");	
	fprintf(stderr,"\n");
	for (ip=0;ip<6;++ip) {
	  for (ic=0;ic<nCoef[i];++ic) {
	    fprintf(stderr,"%1d ",(int)index[i][ic*6 + ip]);
	  }
	  fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
      }
    }

    /*
    ** Calculate all the needed Legendre polynomials up to order lmax
    ** for each of the cosmological input parameters (scaled to [-1,1)).
    */
    //printf("Calculate Legendre Polynomials...\n");
    for (ip=0;ip<6;++ip) {
      Pl[ip] = malloc((lmax+1)*sizeof(double));
      assert(Pl[ip] != NULL);
      x = 2*(CosmoParams[ip] - min[ip])/(max[ip]-min[ip]) - 1.0;
      gsl_sf_legendre_Pl_array(lmax,x,Pl[ip]);
      for (l=0;l<=lmax;++l) Pl[ip][l] *= sqrt(2.0*l + 1.0); // normalize the Pls
      if (bVerbose) {
	for (l=0;l<=lmax;++l) fprintf(stderr,"%1d %.14f\n",l,Pl[ip][l]);
	fprintf(stderr,"\n");
      }
    }

    /*
    ** Now assemble the eigenvalues for the PCAs
    */
    //printf("Assemble EVs for PCAs...\n");
    for (i=0;i<11;++i) {
      lamb[i] = 0;
      for (ic=0;ic<nCoef[i];++ic) {
	/*
	** Below is not very efficient since we uselessly multiply by 1.0
	** most of the time! But as a test it should be OK.
	*/
	prod = 1.0;
	for (ip=0;ip<6;++ip) prod *= Pl[ip][(int)index[i][ic*6 + ip]];
	if (bVerbose && i==5 && coef[i][ic] != 0) {
	    fprintf(stderr,"PCA%02d Product %2d %.14f\n",i,ic,prod);
	    }
	prod *= coef[i][ic];
	lamb[i] += prod;
      }
      if (bVerbose) {
	  fprintf(stderr,"Eigenvalue %1d %.14f\n",i,lamb[i]);
	  }
    }

    /*
    ** Now combine the PCs to form the final log boost (for all z and k).
    ** This is a very data parallel operation. Also if only a specific z 
    ** is required then we should just evaluate it for that redshift.
    ** Some interpolation of redshift and k might also be needed here 
    ** depending on the desired API.
    */

    lnboost = malloc(len_z*nk*sizeof(double));
    boost = malloc(len_z*nk*sizeof(double));

    csmInitialize(&csm);
    csm->val.bComove = 1;
    csm->val.dHubble0 = sqrt(8*M_PI/3);
    csm->val.dOmega0 = CosmoParams[1]/(CosmoParams[3]*CosmoParams[3]);
    csm->val.w0 = CosmoParams[4];
    csm->val.wa = 0.0;      /* this time around we don't consider it */
    csm->val.dOmegaRad = omega_rad/(CosmoParams[3]*CosmoParams[3]);  /* omega_rad (lower case omega) / h^2 */
    csm->val.dOmegaDE = 1.0 - csm->val.dOmega0 - csm->val.dOmegaRad;

    for(zcounter=0; zcounter<len_z; ++zcounter){
   
        z=Redshifts[zcounter];
  
        /* Notice that the interpolation step is only necessary if z != 0 */
        if (z==0) {
            iz = 100;
            for (j=0;j<nk;++j) {
                lnboost[zcounter*nk+j] = mean[iz][j];
            }

            for (i=0;i<11;++i) {
                for (j=0;j<nk;++j) {
                    lnboost[zcounter*nk+j] += lamb[i]*(pca[i + 11*(iz*nk + j)]) ;
                }
            }
        }
        else {
            t0 = csmExp2Time(csm,1.0);
            t200 = csmExp2Time(csm,1.0/201.0); /* time at redshift 200 */
            dDelta = (t0-t200)/nSteps; /* work out the timestep this cosmology would have */
            t = csmExp2Time(csm,1.0/(1.0+z)); /* find the desired proper time in this cosmology */
            dti = (t-t200)/dDelta;  /* find the desired (incl fractional) timestep */
            iz = (int)floor(dti);
            x = dti - iz; /* this is the fractional step */

            //fprintf( stderr, "# requested step for z = %f: %f\n", z,dti);

            for (j=0;j<nk;++j) {
                lnboost[zcounter*nk+j] = ((1-x)*mean[iz][j] + x*mean[iz+1][j]);
            }

            for (i=0;i<11;++i) {
                for (j=0;j<nk;++j) {
                    lnboost[zcounter*nk+j] += lamb[i]*((1-x)*pca[i + 11*(iz*nk + j)] + x*pca[i + 11*((iz+1)*nk + j)]); 
                }
            }
        }
        fprintf(stderr, "# Computed %d of %d boost factors.\n", zcounter+1, len_z);
    } 
       
   for (j=0;j<len_z*nk;++j) {
        boost[j] = exp(lnboost[j]);
        //fprintf(stderr, "boost[%d] = %f\n", j, boost[j]);
    }
      
    *kVals = kvec;
    *Boost = boost;
    *nkVals = nk;
    *nBoost = nk;
    
    int retval;
    retval = munmap (retFid.handle, retFid.size);
    assert(retval == 0);
}
