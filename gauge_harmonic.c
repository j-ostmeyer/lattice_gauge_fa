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

void sample_id(double complex *u, unsigned ns, unsigned nn2, unsigned gd){
	const unsigned dim = ns*nn2, group_dim = gd*gd;
	
	for(unsigned i = 0; i < dim; i++) construct_id(u + i*group_dim, gd);
}

double coupling_fac(double beta, gauge_flags *mode){
	if(mode->gauge_group == 1) return 4 * beta; // factor 2 from definition of generators
	else{
		const unsigned nc = mode->gauge_dim;
		return 2 * beta / nc;
	}
}

double fill_harm_mat(double complex *m, unsigned nn2, unsigned nl, unsigned i, gauge_flags *mode){
	// fill unitary matrix of eigenvectors, return non-zero eigenvalue
	const unsigned loc_dim = nl/2+1, mat_dim = nn2*nn2;
	const double matsubara = M_PI / nl;
	double *w = mode->ddummy;
	double *z = w + nn2;
	double complex *y = m + mat_dim;

	double diag = 0, wmx = 0;
	unsigned ind = i / loc_dim, pos = i % loc_dim, imx = 0;

	for(unsigned k = 0; k < nn2; k++){
		z[k] = matsubara * pos;
		w[k] = sin(z[k]);
		y[k] = cexp(I * z[k]);
		diag += w[k] * w[k];

		if(fabs(w[k]) > wmx) imx = k;

		pos = ind % nl;
		ind /= nl;
	}

	// 0th eigenvector, m*v0 = 0
	const double w_inv = 1./sqrt(diag);
	for(unsigned k = 0; k < nn2; k++) m[k] = w_inv * w[k] * y[k];

	// remaining eigenvectors, m*v_j = diag * v_j
	for(unsigned k = nn2; k < mat_dim; k++) m[k] = 0;
	for(unsigned i = 0, j = 0; i < nn2; i++){
		if(i == imx) continue; // spanning the eigenspace with best pivoting
		j++;
		
		const unsigned shift = j*nn2;
		m[shift + imx] =  y[imx] * w[i];
		m[shift + i  ] = -y[i]   * w[imx];
	}

	// make the matrix unitary
	mod_gram_schmidt(m + nn2, nn2-1, nn2);

	//print_zmat(m, nn2);

	return diag;
}

void mat_mul_basis(double complex *x, double complex *y, double complex *m, unsigned nn2, unsigned ng, unsigned j, int dagger){
	for(unsigned k1 = 0; k1 < nn2; k1++){
		const unsigned shiftM = k1*nn2, shift = k1*ng + j;
		y[shift] = 0;

		for(unsigned k2 = 0; k2 < nn2; k2++){
			if(dagger){
				y[shift] += conj(m[shiftM + k2]) * x[j + k2*ng];
			}else{
				y[shift] += m[k1 + k2*nn2] * x[j + k2*ng];
			}
		}
	}
}

void sample_fourier_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double scale = sqrt(coupling_fac(beta, mode)) / ns, mass2 = mode->fa_mass * mode->fa_mass;
	double complex *m = mode->zdummy;
	double complex *tmp = m + nn2*nn2;

	random_vector(p, dim);

	if(mode->no_fourier_acc) return;

	fftw_execute(fft[0]);

	// project zero eigenmodes to 0
	for(unsigned j = 0; j < nn2*ng; j++) pc[j] *= scale * mode->fa_mass;

	for(unsigned i = 1; i < compl_dim; i++){
		const double ev = fill_harm_mat(m, nn2, nl, i, mode), fac = scale*sqrt(ev + mass2);
		const unsigned shift = i*nn2*ng;

		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			mat_mul_basis(pc + shift, tmp, m, nn2, ng, j, 1);

			tmp[j] *= scale * mode->fa_mass; // 0th eigenvalue is zero
			for(unsigned k = 1; k < nn2; k++) tmp[j + k*ng] *= fac;

			mat_mul_basis(tmp, pc + shift, m, nn2, ng, j, 0);
		}
	}

	fftw_execute(fft[1]);
}

double energy_fourier_momenta(double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double mass2 = mode->fa_mass * mode->fa_mass;
	double complex *m = mode->zdummy;
	double complex *tmp = m + nn2*nn2;
	double en = 0;

	fftw_execute(fft[0]);

	// treat zero eigenmodes
	if(mass2 > 0) for(unsigned j = 0; j < nn2*ng; j++) en += .5/mass2 * norm(pc[j]);

	for(unsigned i = 1; i < compl_dim; i++){
		const unsigned pos = i % loc_dim;
		double weight = .5;

		// elements t=1...Nt/2-1 occur twice (as complex conj pairs), but are stored only once
		if(pos > 0 && pos < (nl+1)/2) weight *= 2;

		const double ev = fill_harm_mat(m, nn2, nl, i, mode), fac = weight/(ev + mass2);
		const unsigned shift = i*nn2*ng;

		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			mat_mul_basis(pc + shift, tmp, m, nn2, ng, j, 1);

			if(mass2 > 0) en += weight/mass2 * norm(tmp[j]); // 0th eigenvalue is zero
			for(unsigned k = 1; k < nn2; k++) en += fac * norm(tmp[j + k*ng]);
		}
	}

	return en / coupling_fac(beta, mode) / ns;
}

void get_fourier_x_dot(double complex *pc, double beta, unsigned ns, unsigned nn2, double h, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim;
	const double scale = h / coupling_fac(beta, mode) / ns, mass2 = mode->fa_mass * mode->fa_mass;
	double complex *m = mode->zdummy;
	double complex *tmp = m + nn2*nn2;

	fftw_execute(fft[0]);

	// project zero eigenmodes to 0 or cutoff 1/m^2
	for(unsigned j = 0; j < nn2*ng; j++) pc[j] *= mass2? scale/mass2 : 0;

	for(unsigned i = 1; i < compl_dim; i++){
		const double ev = fill_harm_mat(m, nn2, nl, i, mode), fac = scale/(ev + mass2);
		const unsigned shift = i*nn2*ng;

		for(unsigned j = 0; j < ng; j++){ // sub-optimal cache-locality, but array should be small enough
			mat_mul_basis(pc + shift, tmp, m, nn2, ng, j, 1);

			tmp[j] *= mass2? scale/mass2 : 0; // potential Gauge-fixing: no/regularised dynamics on zero eigenmodes
			for(unsigned k = 1; k < nn2; k++) tmp[j + k*ng] *= fac;

			mat_mul_basis(tmp, pc + shift, m, nn2, ng, j, 0);
		}
	}

	fftw_execute(fft[2]);
}
