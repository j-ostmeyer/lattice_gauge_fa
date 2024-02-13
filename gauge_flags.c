#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "gauge_flags.h"

void set_flags(FILE *in, gauge_flags *mode){
	mode->gauge_group = 0;
	mode->no_fourier_acc = 0;

	char flag[200];

	while(fscanf(in, "%s\n", flag) != EOF && !mode->gauge_group){
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
		}
	}

	if(!mode->gauge_group) printf("ERROR: Gauge group missing or unknown!\n");
}
