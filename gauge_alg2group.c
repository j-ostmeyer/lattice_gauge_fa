#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifndef MKL_LAPACKE
#include <lapacke.h>
#else
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

#include "gauge_flags.h"
#include "gauge_aux.h"

#include "gauge_alg2group.h"

void alg2group(double *x, double complex *u, double complex *g, double *ev, int dagger, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	double complex *y = g + mat_dim;
	double complex *gr = y + mat_dim;
	double complex *cev = gr + mat_dim;

	//for(unsigned i = 0; i < mat_dim; i++) g[i] = 0;

	switch(mode->gauge_group){
		case 1:
			u[0] = cexp(I * x[0]);
			return;
		case 2:
			get_alg_su2(x, g);
			break;
		case 3:
			get_alg_su3(x, g);
			break;
	}

	for(unsigned j = 0; j < n; j++){
		const unsigned jn = j*n;

		g[jn + j] *= .5; // half diagonal terms

		for(unsigned k = j+1; k < n; k++){
			g[jn + k] *= .5; // half upper triangle terms
			g[k*n + j] = conj(g[jn + k]);
		}
	}

	// diagonalise g
	// nasty hack because zheev doesn't work!
	//LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', n, g, n, ev);
	LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, g, n, cev, NULL, n, gr, n);
	g = gr;
	for(unsigned i = 0; i < n; i++) ev[i] = creal(cev[i]);

	// exp(i ev)
	for(unsigned i = 0; i < n; i++) y[i] = cexp(I * ev[i]);

	// multiply back
	for(unsigned j = 0; j < n; j++){
		const unsigned jn = j*n;

		for(unsigned k = 0; k < n; k++){
			const unsigned kn = k*n;
			double complex tmp = 0;

			for(unsigned i = 0; i < n; i++){
				tmp += g[jn + i] * y[i] * conj(g[kn + i]);
			}

			if(dagger) u[kn + j] = conj(tmp);
			else u[jn + k] = tmp;
		}
	}
}

void get_alg_u1(double *x, double complex *g){
	g[0] = x[0];
}

void get_alg_su2(double *x, double complex *g){
	g[1]  =     x[0]; // sigma_x
	g[1] -= I * x[1]; // sigma_y
	g[0]  =     x[2]; // sigma_z
	g[3]  =    -x[2]; // sigma_z
}

void get_alg_su3(double *x, double complex *g){
	g[1]  =        x[0]; // lambda_1
	g[1] -=    I * x[1]; // lambda_2
	g[0]  =        x[2]; // lambda_3
	g[4]  =       -x[2]; // lambda_3

	g[2]  =        x[3]; // lambda_4
	g[2] -=    I * x[4]; // lambda_5
	g[5]  =        x[5]; // lambda_6
	g[5] -=    I * x[6]; // lambda_7

	g[0] +=   ISQ3*x[7]; // lambda_8
	g[4] +=   ISQ3*x[7]; // lambda_8
	g[8]  =-2*ISQ3*x[7]; // lambda_8
}

void project_tr_lambda(double complex *u, double *x, gauge_flags *mode){
	// trace of u projected onto generators, i.e. x = tr[lambda * u]
	switch(mode->gauge_group){
		case 1:
			project_tr_u1(u, x);
			break;
		case 2:
			project_tr_su2(u, x);
			break;
		case 3:
			project_tr_su3(u, x);
			break;
	}
}

void project_tr_u1(double complex *u, double *x){
	x[0] = creal(u[0]);
}

void project_tr_su2(double complex *u, double *x){
	x[0] = .5 * ( creal(u[1]) + creal(u[2])); // sigma_x
	x[1] = .5 * (-cimag(u[1]) + cimag(u[2])); // sigma_y
	x[2] = .5 * ( creal(u[0]) - creal(u[3])); // sigma_z
}

void project_tr_su3(double complex *u, double *x){
	x[0] = .5 * ( creal(u[1]) + creal(u[3])); // lambda_1
	x[1] = .5 * (-cimag(u[1]) + cimag(u[3])); // lambda_2
	x[2] = .5 * ( creal(u[0]) - creal(u[4])); // lambda_3

	x[3] = .5 * ( creal(u[2]) + creal(u[6])); // lambda_4
	x[4] = .5 * (-cimag(u[2]) + cimag(u[6])); // lambda_5
	x[5] = .5 * ( creal(u[5]) + creal(u[7])); // lambda_6
	x[6] = .5 * (-cimag(u[5]) + cimag(u[7])); // lambda_7

	x[7] = creal(u[0]) + creal(u[4]) - 2*creal(u[8]); // lambda_8
	x[7] *= .5*ISQ3; // lambda_8
}
