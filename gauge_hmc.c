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

void gauge_force(double complex *u, double *p_dot, double beta, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned links = nn/2, ng = mode->num_gen, nc = mode->gauge_dim;
	double complex *pl = mode->zdummy;

	for(unsigned i = 0; i < ns; i++){
		for(unsigned mu = 0; mu < links; mu++){
			const unsigned shift = ng * (mu + i*links);

			sum_of_plaquettes(u, pl, nnt, ns, nn, i, mu, 0, mode);
			project_tr_lambda(pl, p_dot + shift, mode);

			for(unsigned k = 0; k < ng; k++) p_dot[shift + k] *= 2*beta / nc;
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
	const unsigned nd = mode->space_dim, nc = mode->gauge_dim;

	return energyP - beta / nc * pl_av * ns * nd*(nd-1)/2;
}

void update_u(double complex *u, double *p, double complex *pc, double *x_dot, double beta, unsigned ns, unsigned nn, double h, const fftw_plan *fft, gauge_flags *mode){
	const unsigned nn2 = nn/2;
	const unsigned ng = mode->num_gen, dim = ns*nn2;
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	double complex *g = mode->zdummy;
	double complex *up = g + 4*mat_dim;
	double *ev = mode->ddummy;
		
	if(mode->no_fourier_acc){
		const unsigned p_dim = dim*ng;
		for(unsigned k = 0; k < p_dim; k++) x_dot[k] = h*p[k];
	}else get_fourier_x_dot(pc, beta, ns, nn2, h, fft, mode);

	for(unsigned i = 0; i < dim; i++){
		const unsigned shiftP = i*ng, shiftU = i*mat_dim;
		alg2group(x_dot + shiftP, up, g, ev, 0, mode); // exp(i h p)
		mat_mul(up, u + shiftU, g, n); // exp(i h p) * u
		copy_mat(g, u + shiftU, n, 0);
	}
}

void leap_frog(double complex *u, double *p, double complex *pc, double *p_dot, double *x_dot, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, const fftw_plan *fft, gauge_flags *mode){
	const unsigned ng = mode->num_gen, nn2 = nn/2, dim = ns*nn2*ng;
	double h = traj_length/steps;

	update_u(u, p, pc, x_dot, beta, ns, nn, .5*h, fft, mode);

	for(unsigned s = 1; s <= steps; s++){
		gauge_force(u, p_dot, beta, nnt, ns, nn, mode);
		for(unsigned i = 0; i < dim; i++) p[i] += h*p_dot[i];

		if(s == steps) h *= .5; // only half step in the end
		update_u(u, p, pc, x_dot, beta, ns, nn, h, fft, mode);
	}
}

short trajectory(double *p, double complex *pc, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, int iteration, unsigned meas_freq, double *plaquettes_old, const fftw_plan *fft, FILE *res_out, gauge_flags *mode){
	const unsigned nn2 = nn/2, ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim * nn2*ng;
	const unsigned gd = mode->gauge_dim, group_dim = gd*gd, link_dim = ns*nn2*group_dim;
	double *x_dot = p + dim; // dim
	double *p_dot = x_dot + dim; // dim
	double complex *u = pc + compl_dim;
	double complex *u_old = u + link_dim;

	double energy, energyP, plaquettes;
	double energy_old, energyP_old;
	double boltzmann;
	short acc = 1;

	for(unsigned i = 0; i < link_dim; i++) u_old[i] = u[i];

	sample_fourier_momenta(p, pc, beta, ns, nn2, fft, mode);

	energyP_old = energy_momenta(p, pc, beta, ns, nn2, fft, mode);
	energy_old = hamilton(energyP_old, *plaquettes_old, beta, ns, mode);

	leap_frog(u, p, pc, p_dot, x_dot, beta, nnt, ns, nn, steps, traj_length, fft, mode);
	check_unitarity(u, ns, nn, mode);

	energyP = energy_momenta(p, pc, beta, ns, nn2, fft, mode);
	plaquettes = plaquette_av(u, nnt, ns, nn, mode);
	energy = hamilton(energyP, plaquettes, beta, ns, mode);

	//printf("E_p: %g -> %g\n", energyP_old, energyP);
	//printf("plaquette: %g -> %g\n", *plaquettes_old, plaquettes);
	//printf("energies: %g -> %g ,  dE = %g\n", energy_old, energy, energy-energy_old);

	boltzmann = exp(energy_old-energy);
	if(!(boltzmann > genrand64_real2())){ // cumbersome expression to avoid problems with NaN
		// reject (nothing to be done if accepted)
		for(unsigned i = 0; i < link_dim; i++) u[i] = u_old[i];
		plaquettes = *plaquettes_old;
		energyP = energyP_old;
		energy = energy_old;
		acc = 0;
	}

	*plaquettes_old = plaquettes;

	if(iteration % meas_freq == 0){
		double *results = mode->ddummy;
		
		// collection of different results, current dimension = NUM_RES = 7
		results[0] = plaquettes / gd; // average plaquette
		results[1] = plaquettes / gd / plaquettes_old[1] - 1; // relative plaquette deviation from strong coupling result
		results[2] = topo_charge(u, nnt, ns, nn, mode); // topological charge in 2D or 4D
		results[3] = results[2]*results[2]; // squared topological charge
		results[4] = energy; // HMC energy
		results[5] = acc; // acceptance
		results[6] = boltzmann; // exp(-dH), should average to 1

		for(unsigned l = 2; l <= mode->wilson_loops; l++)
			results[NUM_RES+l-2] = wilson_loop_av(u, nnt, ns, nn, l, l, mode) / gd; // average Wilson loops of fixed length

		fprint_results(res_out, results, 1, mode->num_res);
	}

	return acc;
}

void run_hmc(double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, unsigned therm, unsigned meas, unsigned meas_freq, const char *res_name, gauge_flags *mode){

	//printf("Hi from run_hmc with %u flavors!!!\n", flavors); fflush(stdout);
	const unsigned nn2 = nn/2, ng = mode->num_gen, dim = ns*nn2*ng;
	const unsigned nd = mode->space_dim, nl = mode->length_cube, loc_dim = nl/2+1, compl_dim = (ns/nl) * loc_dim * nn2*ng;
	const unsigned gd = mode->gauge_dim, group_dim = gd*gd, link_dim = ns*nn2*group_dim;

	double *p = malloc(3*dim * sizeof(double));
	fftw_complex *pc = (fftw_complex*) fftw_malloc((compl_dim + 2*link_dim) * sizeof(fftw_complex));
	double complex *u = pc + compl_dim;

	mode->zdummy = malloc(20 * (group_dim + ng*nn2*nn2) * sizeof(double complex));
	mode->ddummy = malloc((10 * gd + mode->num_res) * sizeof(double));
	mode->idummy = malloc(nn * sizeof(int));

	sample_id(u, ns, nn2, gd);

	double plaquettes[2];
	plaquettes[0] = plaquette_av(u, nnt, ns, nn, mode);
	plaquettes[1] = strong_coupling_plaquette(beta, mode);

	//printf("%.16g\n", plaquettes[1]);
	//exit(0);

	FILE *res_out = res_name? fopen(res_name, "a"): NULL;

	init_genrand64(time(NULL));

	int *fft_dim = malloc(nd * sizeof(int));
	for(unsigned i = 0; i < nd; i++) fft_dim[i] = nl;
	fftw_plan fft[3]; // p forward/backward, pc -> x_dot
	fft[0] = fftw_plan_many_dft_r2c(nd, fft_dim, nn2*ng, p, NULL, nn2*ng, 1, pc, NULL, nn2*ng, 1, FFTW_MEASURE);
	fft[1] = fftw_plan_many_dft_c2r(nd, fft_dim, nn2*ng, pc, NULL, nn2*ng, 1, p, NULL, nn2*ng, 1, FFTW_MEASURE);
	fft[2] = fftw_plan_many_dft_c2r(nd, fft_dim, nn2*ng, pc, NULL, nn2*ng, 1, p+dim, NULL, nn2*ng, 1, FFTW_MEASURE);
	free(fft_dim);
	
	for(unsigned i = 0; i < therm; i++)
		trajectory(p, pc, beta, nnt, ns, nn, steps, traj_length, i, therm+1, plaquettes, fft, NULL, mode);

	for(unsigned i = 0; i < meas; i++)
		trajectory(p, pc, beta, nnt, ns, nn, steps, traj_length, i, meas_freq, plaquettes, fft, res_out, mode);

	free(p);
	fftw_free(pc);
	free(mode->zdummy);
	free(mode->ddummy);
	free(mode->idummy);
	if(res_out) fclose(res_out);

	for(unsigned i = 0; i < 3; i++)	fftw_destroy_plan(fft[i]);
}
