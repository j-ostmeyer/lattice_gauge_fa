#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "gauge_flags.h"

#include "gauge_lattice.h"

unsigned *construct_lattice(FILE *input, unsigned *ns, unsigned *nn, gauge_flags *mode){
	unsigned *nnt;

	char geometry[100];
	fscanf(input, "geometry = %s\n", geometry);

	if(!strcmp(geometry, "cubic")){
		nnt = construct_cubic(input, ns, nn, mode);
	// TODO Other geometries should be implemented!
	}else{
		printf("The geometry \"%s\" is not known!\n", geometry);
		exit(0);
	}

	return nnt;
}

unsigned *construct_cubic(FILE *input, unsigned *ns, unsigned *nn, gauge_flags *mode){
	unsigned l, dim, d;
	unsigned *N, *nnt;

	fscanf(input, "dimension = %u\n", &dim);
	fscanf(input, "length = %u\n", &l);

	mode->length_cube = l;
	mode->space_dim = dim;

	*nn = 2*dim;
	N = malloc(dim * sizeof(unsigned));
	N[0] = 1;
	for(d = 1; d < dim; d++){
		N[d] = N[d-1]*l;
	}
	*ns = N[d-1]*l;

	const unsigned size = *ns;
	const unsigned neighbours = *nn, links = neighbours/2;
	const unsigned mat_dim = mode->gauge_dim * mode->gauge_dim;

	nnt = malloc(2 * neighbours * size * sizeof(unsigned));

	for(unsigned i = 0; i < size; i++){
		for(d = 0; d < dim; d++){
			const unsigned id = (i/N[d]) % l;

			// nearest sites
			unsigned shift = i*neighbours;
			nnt[shift + d]         = i + ((id+1  )%l - id)*N[d]; // right, top,...
			nnt[shift + d + links] = i + ((id+l-1)%l - id)*N[d]; // left, bottom,...
			const unsigned im = nnt[shift + d + links];

			// nearest links
			shift = (size + i)*neighbours;
			nnt[shift + d]         = mat_dim * (d + links*i); // right, top,...
			nnt[shift + d + links] = mat_dim * (d + links*im); // left, bottom,...
		}
	}

	free(N);
	return nnt;
}
