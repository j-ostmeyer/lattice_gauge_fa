#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include "gauge_flags.h"
#include "gauge_lattice.h"

#include "gauge_aux.h"

double norm(double complex x){
	const double a = creal(x), b = cimag(x);
	return a*a + b*b;
}

void print_vec(double *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\n", vec[i]);
}

void print_mat(double *m, unsigned long n){
	unsigned i, k;

	printf("{");
	for(i = 0; i < n; i++){
		printf("{");
		for(k = 0; k < n; k++, m++) printf("%.15g,\t", *m);
		printf("},\n");
	}
	printf("}\n");
}

void fprint_results(FILE *out, double *res, unsigned n, unsigned r){
	if(out){
		for(unsigned i = 0; i < n; i++){
			for(unsigned k = 0; k < r; k++, res++) fprintf(out, "%.15g\t", *res);
			fprintf(out, "\n");
		}
		fflush(out);
	}
}

double average(double *x, unsigned n){
	double res = 0;

	for(unsigned i = 0; i < n; i++) res += x[i];

	return res/n;
}

double average_sq(double *x, unsigned n){
	double res = 0;

	for(unsigned i = 0; i < n; i++) res += pow(x[i], 2);

	return res/n;
}

double scalar_dot(double *x, double *y, unsigned d){
	double total = 0;

	for(unsigned i = 0; i < d; i++) total += x[i]*y[i];

	return total;
}

double reTr(double complex *u, unsigned d){
	const unsigned d1 = d + 1;
	double tr = 0;

	for(unsigned i = 0; i < d; i++) tr += creal(u[i*d1]);

	return tr;
}

void construct_id(double complex *m, unsigned n){
	// construct identity
	const unsigned n2 = n*n;
	const unsigned n1 = n+1;

	for(unsigned i = 0; i < n2; i++) m[i] = 0;
	for(unsigned i = 0; i < n; i++) m[i*n1] = 1;
}

void dagger(double complex *m, unsigned n){
	for(unsigned i = 0; i < n; i++){
		const unsigned shift = i*n;

		for(unsigned j = i+1; j < n; j++){
			const unsigned pos = shift + j;
			const unsigned posT = j*n + i;

			const double complex tmp = m[pos];
			m[pos] = conj(m[posT]);
			m[posT] = conj(tmp);
		}
	}
}

void dagger_sym(double complex *m, unsigned n){
	for(unsigned i = 0; i < n; i++){
		const unsigned shift = i*n;

		for(unsigned j = i+1; j < n; j++){
			const unsigned pos = shift + j;
			const unsigned posT = j*n + i;

			const double complex tmp = .5*(m[pos] + conj(m[posT]));
			m[pos] = tmp;
			m[posT] = tmp;
		}
	}
}

void dagger_asym(double complex *m, unsigned n){
	for(unsigned i = 0; i < n; i++){
		const unsigned shift = i*n;

		for(unsigned j = i+1; j < n; j++){
			const unsigned pos = shift + j;
			const unsigned posT = j*n + i;

			const double complex tmp = .5*I*(m[pos] - conj(m[posT]));
			m[pos] = tmp;
			m[posT] = tmp;
		}
	}
}

void mat_mul(double complex *x, double complex *y, double complex *res, unsigned n){
	for(unsigned j = 0; j < n; j++){
		const unsigned jn = j*n;

		for(unsigned k = 0; k < n; k++) res[jn + k] = 0;

		for(unsigned i = 0; i < n; i++){
			const unsigned in = i*n;
			const double complex tmp = x[jn + i];

			for(unsigned k = 0; k < n; k++){
				res[jn + k] += tmp * y[in + k];
			}
		}
	}
}

void mat_prod(double complex *u, unsigned n, unsigned d){
	// product of n d-dimensional matrices, e.g. for staple
	// right-most matrix first in array u, working space last
	
	const unsigned mat_dim = d*d;
	double complex *prod1 = u + (n + n%2)*mat_dim;
	double complex *prod2 = u + (n + 1 - n%2)*mat_dim;

	switch(n){
		case 0:
			construct_id(prod1, d);
			return;
		case 1:
			for(unsigned i = 0; i < mat_dim; i++) prod2[i] = u[i];
			return;
		default:
			mat_mul(u + mat_dim, u, prod1, d);
	}

	for(unsigned k = 2; k < n; k++){
		mat_mul(u + k*mat_dim, prod1, prod2, d);

		double complex *dummy = prod1;
		prod1 = prod2;
		prod2 = dummy;
	}
}

double trace_prod(double complex *u, unsigned n, unsigned d){
	// trace of product of n d-dimensional matrices, e.g. for plaquette
	// right-most matrix first in array u, working space last
	
	const unsigned mat_dim = d*d;
	double complex *prod;

	switch(n){
		case 0:
			return d;
		case 1:
			return reTr(u, d);
		case 2:
			prod = u + mat_dim;
			break;
		default:
			prod = u + n*mat_dim;
			mat_prod(u + mat_dim, n-1, d);
	}

	double tr = 0;

	for(unsigned k = 0; k < d; k++){
		const unsigned kd = k*d;
		for(unsigned i = 0; i < d; k++){
			const double complex a = prod[kd + i], b = u[i*d + k];
			tr += creal(a) * creal(b) - cimag(a) * cimag(b);
		}
	}

	return tr;
}