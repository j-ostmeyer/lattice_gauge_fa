#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "gauge_flags.h"
#include "gauge_aux.h"
#include "gauge_lattice.h"

#include "gauge_hmc.h"

int main(int argc, char **argv){
	switch(argc){
		case 3:
			break;
		default:
			printf("Error: One input and one output file needed!\n");
			return 0;
	}

	FILE *in = fopen(argv[1], "r");

	double beta;
	if(fscanf(in, "beta = %lg\n", &beta) != 1){
		printf("Error: Need beta!\n");
		return 0;
	}

	unsigned steps;
	double traj_length;
	if(fscanf(in, "MD steps, traj. length = %u, %lg\n", &steps, &traj_length) != 2){
		printf("Error: Need steps, traj-length!\n");
		return 0;
	}

	unsigned therm, meas, meas_freq;
	if(fscanf(in, "therm. steps = %u\nmeas. steps, freq. = %u, %u\n", &therm, &meas, &meas_freq) != 3){
		printf("Error: Need therm, meas, meas-freq!\n");
		return 0;
	}

	//printf("%g, %u, %g, %u, %u, %u\n", beta, steps, traj_length, therm, meas, meas_freq);

	gauge_flags mode[1];
	set_flags(in, mode);

	unsigned ns, nn;
	unsigned *nnt = construct_lattice(in, &ns, &nn, mode);

	fclose(in);

	if(traj_length == 0){
		traj_length = HALF_M_PI; // de-correlate all phonon modes exactly with FA
	}

	run_hmc(beta, nnt, ns, nn, steps, traj_length, therm, meas, meas_freq, argv[2], mode);

	free(nnt);

	return 0;
}
