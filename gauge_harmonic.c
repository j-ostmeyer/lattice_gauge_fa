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

#include "mt19937-64.h"
#include "gauge_random.h"

#include "gauge_flags.h"
#include "gauge_aux.h"

#include "gauge_harmonic.h"

void fill_harm_mat(double complex *u, double beta, unsigned nn2, unsigned nl, unsigned i, gauge_flags *mode){
	// TODO this only works for hyper-cubic lattices!
	const unsigned loc_dim = nl/2+1;
	const double matsubara = M_PI / nl;
	double *w = mode->ddummy;
	double *z = w + nn2;

	double diag = 0;
	unsigned ind = i / loc_dim, pos = i % loc_dim;

	for(unsigned k = 0; k < nn2; k++){
		z[k] = matsubara * pos;
		w[k] = sin(z[k]);
		diag += w[k] * w[k];
		pos = ind % nl;
		ind /= nl;
	}

	diag *= beta; // factor ns because fft accumulates it from there-back trafo
	for(unsigned k = 0; k < nn2; k++) u[k*(nn2+1)] = diag;

	for(unsigned k1 = 1; k1 < nn2; k1++){
		const unsigned shift = k1*nn2;
		for(unsigned k2 = 0; k2 < k1; k2++){
			u[shift + k2] = beta * w[k1] * w[k2] * cexp(I * (z[k1] - z[k2]));
		}
	}
}

void sample_harmonic(double *x, double complex *xc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	// TODO implement properly
	const unsigned ng = mode->num_gen, dim = ns*nn2*ng;
	
	for(unsigned i = 0; i < dim; i++) x[i] = 0;
}

void sample_fourier_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	double complex *u = mode->zdummy;

	random_vector(p, dim);

	if(mode->no_fourier_acc) return;

	fftw_execute(fft[1]);

	const double pii = 1. / M_PI / ns;
	for(unsigned j = 0; j < nn2*ng; j++) pc[j] *= pii;

	for(unsigned i = 1; i < compl_dim; i++){
		fill_harm_mat(u, beta*ns, nn2, nl, i, mode);
		LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'L', nn2, u, nn2); // Cholesky-decompose

		const unsigned shift = i*nn2*ng;
		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			const unsigned shiftP = shift + j;
			
			for(unsigned k1 = 0; k1 < nn2; k1++){
				const unsigned shiftU = k1*nn2;
				double complex tmp = 0;

				for(unsigned k2 = 0; k2 <= k1; k2++){
					tmp += u[shiftU + k2] * pc[shiftP + k2*ng];
				}
				pc[shiftP + k1*ng] = tmp;
			}
		}
	}

	fftw_execute(fft[3]);
}

double energy_fourier_momenta(double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	double complex *u = mode->zdummy;
	double en = 0;

	fftw_execute(fft[1]);

	const double pi2 = M_PI*M_PI / ns;
	for(unsigned j = 0; j < nn2*ng; j++) en += norm(pc[j]) * pi2;

	for(unsigned i = 1; i < compl_dim; i++){
		fill_harm_mat(u, beta*ns, nn2, nl, i, mode);
		LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'L', nn2, u, nn2); // Cholesky-decompose
		LAPACKE_zpotri(LAPACK_ROW_MAJOR, 'L', nn2, u, nn2); // invert

		const unsigned shift = i*nn2*ng;
		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			const unsigned shiftP = shift + j;
			
			for(unsigned k1 = 0; k1 < nn2; k1++){
				const unsigned shiftU = k1*nn2;
				double complex tmp = 0;

				unsigned k2;
				for(k2 = 0; k2 <= k1; k2++)
					tmp += u[shiftU + k2] * pc[shiftP + k2*ng];
				for(; k2 < nn2; k2++)
					tmp += conj(u[k2*nn2 + k1]) * pc[shiftP + k2*ng];

				const double complex proj = pc[shiftP + k1*ng];
				en += creal(proj)*creal(tmp) + cimag(proj)*cimag(tmp);
			}
		}
	}

	return en;
}

void evolve_fields(double complex *xc, double beta, unsigned ns, unsigned nn2, double h, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double ivol = 1. / ns;
	double complex *u = mode->zdummy;
	double complex *pc = xc + dim;

	fftw_execute(fft[0]);
	fftw_execute(fft[1]);

	const double pi2 = M_PI*M_PI * h;
	for(unsigned j = 0; j < nn2*ng; j++){
		xc[j] += pi2 * pc[j];
		xc[j] *= ivol;
	}

	for(unsigned i = 1; i < compl_dim; i++){
		fill_harm_mat(u, beta*ns, nn2, nl, i, mode);
		LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'L', nn2, u, nn2); // Cholesky-decompose
		LAPACKE_zpotri(LAPACK_ROW_MAJOR, 'L', nn2, u, nn2); // invert

		const unsigned shift = i*nn2*ng;
		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			const unsigned shiftP = shift + j;
			
			for(unsigned k1 = 0; k1 < nn2; k1++){
				const unsigned shiftU = k1*nn2;
				double complex tmp = 0;

				unsigned k2;
				for(k2 = 0; k2 <= k1; k2++)
					tmp += u[shiftU + k2] * pc[shiftP + k2*ng];
				for(; k2 < nn2; k2++)
					tmp += conj(u[k2*nn2 + k1]) * pc[shiftP + k2*ng];

				const unsigned pos = shiftP + k1*ng;
				xc[pos] += h * pc[pos];
				xc[pos] *= ivol;
			}
		}
	}

	fftw_execute(fft[2]);
}
