#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>

#include "mt19937-64.h"
#include "gauge_random.h"

#include "gauge_flags.h"
#include "gauge_aux.h"
#include "gauge_lattice.h"
#include "gauge_harmonic.h"
#include "gauge_alg2group.h"
#include "gauge_plaquette.h"

#include "gauge_hmc.h"

void gauge_force(double *x, double *p_dot, double beta, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2, ng = mode->num_gen;
	double complex *pl = mode->zdummy;
	double complex *g = pl + mat_dim;
	double complex *u = g + mat_dim;
	double *ev = mode->ddummy;

	for(unsigned i = 0; i < ns; i++){
		for(unsigned mu = 0; mu < links; mu++){
			const unsigned shift = ng * (mu + i*links);

			sum_of_plaquettes(x, pl, u, g, ev, nnt, ns, nn, i, mu, 0, mode);
			project_tr_lambda(pl, p_dot + shift, mode);

			for(unsigned k = 0; k < ng; k++) p_dot[shift + k] *= beta;
		}
	}
}

double energy_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode){
	if(mode->no_fourier_acc){
		const unsigned ng = mode->num_gen, dim = ns*nn2*ng;
		return .5 * dim * average_sq(p, dim);
	}else{
		return energy_fourier_momenta(pc, beta, ns, nn2, fft, mode);
	}
}

double hamilton(double energyP, double pl_av, double beta, unsigned ns, gauge_flags *mode){
	const unsigned nd = mode->space_dim;

	return energyP - .5 * beta * pl_av * ns * nd*(nd-1)/2;
}

void update_x(double *x, double *p, double complex *xc, double beta, unsigned ns, unsigned nn, double h, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, nn2 = nn/2, dim = ns*nn2*ng;

	if(mode->no_fourier_acc){
		for(unsigned i = 0; i < dim; i++) x[i] += h*p[i];
	}else evolve_fields(xc, beta, ns, nn2, h, fft, mode);
}

void leap_frog(double *x, double *p, double *p_dot, double complex *xc, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, nn2 = nn/2, dim = ns*nn2*ng;
	double h = traj_length/steps;

	update_x(x, p, xc, beta, ns, nn, .5*h, fft, mode);

	for(unsigned s = 1; s <= steps; s++){
		gauge_force(x, p_dot, beta, nnt, ns, nn, mode);
		for(unsigned i = 0; i < dim; i++) p[i] += h*p_dot[i];

		if(s == steps) h *= .5; // only half step in the end
		update_x(x, p, xc, beta, ns, nn, h, fft, mode);
	}
}

short trajectory(double *x, double complex *xc, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, int iteration, unsigned meas_freq, double *plaquettes_old, const fftw_plan *fft, FILE *res_out, gauge_flags *mode){
	const unsigned nn2 = nn/2, ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim * nn2*ng;
	double *p = x + dim; // dim
	double *p_dot = p + dim; // dim
	double *x_old = p_dot + dim; // dim
	double complex *pc = xc + compl_dim;

	double energy, energyP, plaquettes;
	double energy_old, energyP_old;
	double boltzmann;
	short acc = 1;

	for(unsigned i = 0; i < dim; i++) x_old[i] = x[i];

	sample_fourier_momenta(p, pc, beta, ns, nn2, fft, mode);

	energyP_old = energy_momenta(p, pc, beta, ns, nn2, fft, mode);
	energy_old = hamilton(energyP_old, *plaquettes_old, beta, ns, mode);

	leap_frog(x, p, p_dot, xc, beta, nnt, ns, nn, steps, traj_length, fft, mode);

	energyP = energy_momenta(p, pc, beta, ns, nn2, fft, mode);
	plaquettes = plaquette_av(x, nnt, ns, nn, mode);
	energy = hamilton(energyP, plaquettes, beta, ns, mode);

	//printf("logDetM: %g -> %g\n", logDetM_old, logDetM);
	//printf("energies: %g -> %g ,  dE = %g\n", energy_old, energy, energy-energy_old);

	boltzmann = exp(energy_old-energy);
	if(!(boltzmann > genrand64_real2())){ // cumbersome expression to avoid problems with NaN
		// reject (nothing to be done if accepted)
		for(unsigned i = 0; i < dim; i++) x[i] = x_old[i];
		plaquettes = *plaquettes_old;
		energyP = energyP_old;
		energy = energy_old;
		acc = 0;
	}

	*plaquettes_old = plaquettes;

	if(iteration % meas_freq == 0){
		double *results = mode->ddummy;
		
		// collection of different results, current dimension = NUM_RES = 5
		results[0] = plaquettes; // average plaquette
		results[1] = plaquettes - strong_coupling_plaquette(beta, mode); // plaquette deviation from strong coupling result
		results[2] = energy; // HMC energy
		results[3] = acc; // acceptance
		results[4] = boltzmann; // exp(-dH), should average to 1

		fprint_results(res_out, results, 1, NUM_RES);
	}

	return acc;
}

void run_hmc(double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, unsigned therm, unsigned meas, unsigned meas_freq, const char *res_name, gauge_flags *mode){

	//printf("Hi from run_hmc with %u flavors!!!\n", flavors); fflush(stdout);
	const unsigned nn2 = nn/2, ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nd = mode->space_dim, nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim * nn2*ng;
	const unsigned gd = mode->gauge_dim, group_dim = gd*gd;

	double *x = malloc(4*dim * sizeof(double));
	fftw_complex *xc = (fftw_complex*) fftw_malloc(2*compl_dim * sizeof(fftw_complex));

	mode->zdummy = malloc(10 * (group_dim + nn2*nn2) * sizeof(double complex));
	mode->ddummy = malloc(10 * gd * sizeof(double));
	mode->idummy = malloc(nn * sizeof(int));

	FILE *res_out = res_name? fopen(res_name, "a"): NULL;

	init_genrand64(time(NULL));

	int *fft_dim = malloc(nd * sizeof(int));
	for(unsigned i = 0; i < nd; i++) fft_dim[i] = nl;
	fftw_plan fft[4]; // x/p forward/backward in every combination
	fft[0] = fftw_plan_many_dft_r2c(nd, fft_dim, nn2*ng, x, NULL, nn2*ng, 1, xc, NULL, nn2*ng, 1, FFTW_MEASURE);
	fft[1] = fftw_plan_many_dft_r2c(nd, fft_dim, nn2*ng, x+dim, NULL, nn2*ng, 1, xc+compl_dim, NULL, nn2*ng, 1, FFTW_MEASURE);
	fft[2] = fftw_plan_many_dft_c2r(nd, fft_dim, nn2*ng, xc, NULL, nn2*ng, 1, x, NULL, nn2*ng, 1, FFTW_MEASURE);
	fft[3] = fftw_plan_many_dft_c2r(nd, fft_dim, nn2*ng, xc+compl_dim, NULL, nn2*ng, 1, x+dim, NULL, nn2*ng, 1, FFTW_MEASURE);
	free(fft_dim);

	sample_harmonic(x, xc, beta, ns, nn2, fft, mode);

	double plaquettes = plaquette_av(x, nnt, ns, nn, mode);
	
	for(unsigned i = 0; i < therm; i++)
		trajectory(x, xc, beta, nnt, ns, nn, steps, traj_length, i, therm+1, &plaquettes, fft, NULL, mode);

	for(unsigned i = 0; i < meas; i++)
		trajectory(x, xc, beta, nnt, ns, nn, steps, traj_length, i, meas_freq, &plaquettes, fft, res_out, mode);

	free(x);
	fftw_free(xc);
	free(mode->zdummy);
	free(mode->ddummy);
	free(mode->idummy);
	if(res_out) fclose(res_out);

	for(unsigned i = 0; i < 4; i++)	fftw_destroy_plan(fft[i]);
}
