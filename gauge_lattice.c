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
	// Geometries below are not adjusted and will fail!
	//}else if(!strcmp(geometry, "rectangular")){
	//	nnt = construct_rectangular(input, ns, nn, s);
	//}else if(!strcmp(geometry, "triangular")){
	//	nnt = construct_triangular(input, ns, nn, s);
	//}else if(!strcmp(geometry, "alltoall")){
	//	nnt = construct_alltoall(input, ns, nn, s);
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
	const unsigned generators = mode->num_gen;

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
			nnt[shift + d]         = generators * (d + links*i); // right, top,...
			nnt[shift + d + links] = generators * (d + links*im); // left, bottom,...
		}
	}

	free(N);
	return nnt;
}

unsigned *construct_rectangular(FILE *input, unsigned *ns, unsigned *nn, double **s){
	unsigned dim, d, i;
	unsigned *l, *N, *nnt;
	char boundaries[100];
	int id;

	fscanf(input, "dimension=%u\n", &dim);
	l = malloc(dim * sizeof(unsigned));
	N = malloc(dim * sizeof(unsigned));
	fscanf(input, "lengths=%u ", l);
	for(d = 1; d < dim; d++){
		fscanf(input, "%u ", l+d);
	}
	fscanf(input, "boundaries=%s\n", boundaries);

	*nn = 2*dim;
	N[0] = 1;
	for(d = 1; d < dim; d++){
		N[d] = N[d-1]*l[d-1];
	}
	*ns = N[d-1]*l[d-1];

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc((size+1) * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	if(!strcmp(boundaries, "periodic")){
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l[d];
				nnt[i*neighbours + 2*d] = i + ((id+1)%l[d] - id)*N[d];
				nnt[i*neighbours + 2*d + 1] = i + ((id+l[d]-1)%l[d] - id)*N[d];
			}
		}
	}else if(!strcmp(boundaries, "open") || !strcmp(boundaries, "closed")){
		if(!strcmp(boundaries, "open")) (*s)[size] = 0;
		else (*s)[size] = 1;
		for(i = 0; i < size; i++){
			for(d = 0; d < dim; d++){
				id = (i/N[d]) % l[d];
				if(id < l[d]-1) nnt[i*neighbours + 2*d] = i + N[d];
				else nnt[i*neighbours + 2*d] = size;
				if(id > 0) nnt[i*neighbours + 2*d + 1] = i - N[d];
				else nnt[i*neighbours + 2*d + 1] = size;
			}
		}
	}else{
		printf("Boundary conditions \"%s\" not known!\n", boundaries);
		exit(0);
	}

	free(l);
	free(N);
	return nnt;
}

unsigned *construct_triangular(FILE *input, unsigned *ns, unsigned *nn, double **s){
	// constructs a parallelogram with periodic boundaries
	unsigned l1, l2, i;
	unsigned *nnt;
	int i1, i2;

	fscanf(input, "lengths=%u %u\n", &l1, &l2);

	*nn = 6;
	*ns = l1*l2;

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc(size * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	for(i = 0; i < size; i++){
		i1 = i%l1;
		i2 = i/l1;
		nnt[i*neighbours] = i + (i1+1)%l1 - i1; // right
		nnt[i*neighbours+1] = i + (i1+l1-1)%l1 - i1; // left
		nnt[i*neighbours+2] = i + ((i2+1)%l2 - i2)*l1; // top right
		nnt[i*neighbours+3] = i + ((i2+l2-1)%l2 - i2)*l1; // bottom left
		nnt[i*neighbours+4] = nnt[i*neighbours+2] + (i1+l1-1)%l1 - i1; // top left
		nnt[i*neighbours+5] = nnt[i*neighbours+3] + (i1+1)%l1 - i1; // bottom right
	}

	return nnt;
}

unsigned *construct_alltoall(FILE *input, unsigned *ns, unsigned *nn, double **s){
	// connects every point to every other point
	unsigned i, pos, k;
	unsigned *nnt;

	fscanf(input, "size=%u\n", ns);

	*nn = *ns-1;

	const unsigned size = *ns;
	const unsigned neighbours = *nn;
	*s = malloc(size * sizeof(double));
	nnt = malloc(neighbours * size * sizeof(unsigned));

	for(i = 0, pos = 0; i < size; i++){
		for(k = 0; k < i; k++, pos++) nnt[pos] = k;
		for(k = i+1; k < size; k++, pos++) nnt[pos] = k;
	}

	return nnt;
}
