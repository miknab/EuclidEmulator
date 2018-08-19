/* This file is part of EuclidEmulator 
 * Copyright (c) 2018 Mischa Knabenhans
 *
 * EuclidEmulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EuclidEmulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <EuclidEmulator.h>

int main(int argc, char* argv[]){

    fprintf(stderr, "# EuclidEmulator Version %d.%d Copyright (C) 2018 Mischa Knabenhans & Joachim Stadel\n",
            EUCLID_EMULATOR_VERSION_MAJOR,
            EUCLID_EMULATOR_VERSION_MINOR);
    fprintf(stderr, "# This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.\n");
    fprintf(stderr, "# This is free software, and you are welcome to redistribute it\n");
    fprintf(stderr, "# under certain conditions; type `show c' for details.\n");

    // Check number of input arguments: 8 arguments are expected:
    // 1 executable name + 6 cosmological parameter + 1 redshift value
    assert(argc>=8);

    // Convert input arguments into appropriate format
    double omega_b = atof(argv[1]);
    double omega_m = atof(argv[2]);
    double n_s = atof(argv[3]);
    double h = atof(argv[4]);
    double w_0 = atof(argv[5]);
    double sigma_8 = atof(argv[6]);
    double* zvec;

    double cosmo_params[6] = {omega_b, omega_m, n_s, h, w_0, sigma_8};
    
    int nz = argc-7; 
    zvec = malloc( nz*sizeof(double) );
    int counter;
    for(counter = 0; counter < nz; counter++){
        zvec[counter] = atof(argv[7+counter]);
        assert(zvec[counter]<=5.0); 
    }
 
    // The radiation density is fixed by the CMB temperature (which is 
    // fixed as well for our purposes).
    const double omega_rad = 4.183709411969527e-5; /* corresponds to 2.755 K Tcmb */


    fprintf( stderr, "#\n# Cosmology:\n");
    fprintf( stderr, "#\tdOmega0: %.8f\n", omega_m/(h*h));
    fprintf( stderr, "#\tdOmegaRad: %.8f\n", omega_rad/(h*h));
    fprintf( stderr, "#\tdOmegaDE: %.8f\n\n", 1.-omega_m/(h*h)-omega_rad/(h*h));

    // Prepare data containers for return values of emulator function
    double* k_values;      // will contain the k modes at which the 
                            // emulator is evaluated
    int length_k_values;    // will contain the number of k modes

    double* boost;         // will contain the boost factor values 
                            // for the k modes in k_values
    int length_boost;       // will contain the number of boost values

    /* REMARK 1:
    ** We declare these data containers and pass them to the EucEmu function
    ** This allows us to be ignorant about the length of the arrays returned
    ** by EucEmu. We pass a pointer to a pointer to the first element of the
    ** two arrays k_values and boost. This allows for modification of the 
    ** pointers pointing directly to the array entries (because the pointer
    ** passed to and returned by the function EucEmu is not altered). Inside
    ** EucEmu the intermediate pointer iterates over all array cells filling
    ** the array with the corresponding values.
    ** Similarly, as we assign a values to length_k_values and length_boost
    ** only inside EucEmu, we modify the content of these two variables inside
    ** the function. This is only allowed by passing by reference which is 
    ** exactly what is done.
    */

    /* REMARK 2: 
    ** length_k_values and length_boost are supposed to be equal. We return
    ** them from the EucEmu function for sanity checking
    */

    // Call the emulator
    EucEmu(&cosmo_params[0], &zvec[0], nz, &k_values, &length_k_values, &boost, &length_boost);

    assert(length_k_values == length_boost);

    // Print results to stdout
    fprintf(stderr,"# ===============================\n");
    printf("#k\t");
    int zcounter;
    for (zcounter=0;zcounter<nz;zcounter++){
	printf("\tz=%.4f", zvec[zcounter]);
    }
    printf("\n#\n");
    int j;
    for(j=0; j<length_k_values; j++){
	printf("%f", k_values[j]);
        for(zcounter=0;zcounter<nz;zcounter++){
	    printf("\t%f", boost[zcounter*length_k_values+j]);
    	}
	printf("\n");
    }

    fprintf(stderr, "Done!\n");

    return 0;
}
