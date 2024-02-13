#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "mt19937-64.h"
#include "gauge_random.h"

#include "gauge_aux.h"

void random_vector(double *x, unsigned N){
	// Box-Muller algorithm for normal distribution
	// quick and dirty implementation, not parallelisable!

	for(unsigned k = 0; k < N; k += 2){
		const double u1 = sqrt(-2 * log(genrand64_real2()));
		const double u2 = M_2PI * genrand64_real2();

		x[k] = u1 * cos(u2);
		x[(k+1)%N] = u1 * sin(u2); // modulo only relevant if N odd
	}
}

void random_cvector(double complex *x, unsigned N){
	// Box-Muller algorithm for normal distribution
	// quick and dirty implementation, not parallelisable!

	for(unsigned k = 0; k < N-1; k++){
		const double u1 = sqrt(-2 * log(genrand64_real2()));
		const double u2 = M_2PI * genrand64_real2();

		if(k){
			x[k] = u1 * (cos(u2) + I*sin(u2));
		}else{
			x[0] = u1 * cos(u2); // first element has to be real
			x[N-1] = u1 * sin(u2); // last element has to be real
		}
	}
}

