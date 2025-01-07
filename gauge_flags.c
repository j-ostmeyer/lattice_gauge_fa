#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "gauge_flags.h"

#include "gauge_aux.h"

void set_flags(FILE *in, gauge_flags *mode){
	mode->gauge_group = 0;
	mode->no_fourier_acc = 0;
	mode->wilson_loops = 1;
	mode->fa_mass = 0;

	char flag[200];

	while(!mode->gauge_group){
		fscanf(in, "%s\n", flag);
		//printf("%s\n", flag);

		if(!strcmp(flag, "U1")){
			mode->gauge_group = 1;
			mode->gauge_dim = 1;
			mode->num_gen = 1;
		}else if(!strcmp(flag, "SU2")){
			mode->gauge_group = 2;
			mode->gauge_dim = 2;
			mode->num_gen = 3;
		}else if(!strcmp(flag, "SU3")){
			mode->gauge_group = 3;
			mode->gauge_dim = 3;
			mode->num_gen = 8;
		}else if(!strcmp(flag, "NO_FOURIER_ACC")){
			mode->no_fourier_acc = 1;
		}else if(!strcmp(flag, "FA_MASS")){
			fscanf(in, "%lg\n", &mode->fa_mass);
		}else if(!strcmp(flag, "WILSON_LOOPS")){
			fscanf(in, "%u\n", &mode->wilson_loops);
		}
	}

	mode->num_res = NUM_RES + mode->wilson_loops-1;

	if(!mode->gauge_group) printf("ERROR: Gauge group missing or unknown!\n");
}
