#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/io.h>
#include <sys/mman.h>
#include <assert.h>
#include <gsl/gsl_sf_legendre.h>
#include "cosmo.h"

// Struct for return type of EucEmu()
struct FID{
    double *handle;
    int size;
};

// Actual emulator function
//void EucEmu(double *CosmoParams, double z);
void EucEmu(double *CosmoParams, double z, double **kVals, int *nkVals, double **Boost, int *nBoost);
