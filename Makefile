all: EuclidEmulator.c cosmo.c cosmo.h
	cc -O3 -g -L. -o ee EuclidEmulator.c cosmo.c -lgsl -lgslcblas -lm

clean:
	rm ee
