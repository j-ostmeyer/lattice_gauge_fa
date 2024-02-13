#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "gauge_flags.h"
#include "gauge_aux.h"
#include "gauge_alg2group.h"

#include "gauge_plaquette.h"

double plaquette_av(double *x, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned nd = mode->space_dim;
	double complex *g = mode->zdummy;
	double complex *u = g + 4*mat_dim;
	double *ev = mode->ddummy;
	double plaquettes = 0;

	for(unsigned i = 0; i < ns; i++){
		for(unsigned nu = 1; nu < links; nu++){
			for(unsigned mu = 0; mu < nu; mu++){
				plaquettes += plaquette_fields(x, u, g, ev, nnt, ns, nn, i, mu, nu, mode);
			}
		}
	}

	return plaquettes * 2/nd/(nd-1) / ns;
}

double plaquette_fields(double *x, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	const unsigned *nnl = nnt + ns*nn;

	alg2group(x + nnl[pos*nn + mu], u + 3*mat_dim, g, ev, 0, mode);
	alg2group(x + nnl[nnt[pos*nn + mu]*nn + nu], u + 2*mat_dim, g, ev, 0, mode);
	alg2group(x + nnl[nnt[pos*nn + nu]*nn + mu], u + 1*mat_dim, g, ev, 1, mode);
	alg2group(x + nnl[pos*nn + nu], u, g, ev, 1, mode);

	return trace_prod(u, 4, n);
}

void staple_fields(double *x, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned *nnl = nnt + ns*nn;

	const int turn = nu >= links;

	alg2group(x + nnl[nnt[pos*nn + mu]*nn + nu], u + 2*mat_dim, g, ev, turn, mode);
	alg2group(x + nnl[nnt[pos*nn + nu]*nn + mu], u + 1*mat_dim, g, ev, 1, mode);
	alg2group(x + nnl[pos*nn + nu], u, g, ev, 1 - turn, mode);

	mat_prod(u, 3, n);
}

void sum_of_plaquettes(double *x, double complex *pl, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, int direction, gauge_flags *mode){
	// calculates the sum of all plaquettes around U_mu(n)
	// direction:  1 only for U_{mu,nu}(n)
	// direction: -1 only for U_{mu,nu}(n)^+
	// direction:  0 all for i/2*(U_{mu,nu}(n) - U_{mu,nu}(n)^+)
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	double complex *st = u + 6*mat_dim;
	const unsigned *nnl = nnt + ns*nn;

	sum_of_staples(x, st, u, g, ev, nnt, ns, nn, pos, mu, mode);
	alg2group(x + nnl[pos*nn + mu], u + 3*mat_dim, g, ev, 0, mode);

	mat_mul(u + 3*mat_dim, st, pl, n);

	if(direction == -1) dagger(pl, n);
	if(direction == 0) dagger_asym(pl, n);
}

void sum_of_staples(double *x, double complex *st, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, gauge_flags *mode){
	// calculates the sum of all staples around U_mu(n)
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	double complex *prod = u + 3*mat_dim;

	for(unsigned i = 0; i < mat_dim; i++) st[i] = 0;
	
	for(unsigned nu = 0; nu < nn; nu++){
		if(nu % links == mu % links) continue;

		staple_fields(x, u, g, ev, nnt, ns, nn, pos, mu, nu, mode);

		for(unsigned i = 0; i < mat_dim; i++) st[i] += prod[i];
	}
}

double strong_coupling_plaquette(double beta, gauge_flags *mode){
	const unsigned d = mode->space_dim;
	switch(mode->gauge_group){
		case 1:
			return 2*strong_coupling_pl_u1(beta, d);
		case 2:
			return 2*strong_coupling_pl_su2(beta, d);
		case 3:
			return 2*strong_coupling_pl_su3(beta, d);
		default:
			return 0;
	}
}

double strong_coupling_pl_u1(double b, unsigned d){
	double p = b / 8.;
	p -= pow(b, 3) / 256.;
	p += pow(b, 5) * (d / 1024. - 11 / 6144.);
	p -= pow(b, 7) * (d / 4096. - 757 / 1572864.);

	return p;
}

double strong_coupling_pl_su2(double b, unsigned d){
	double p = b / 4.;
	p -= pow(b, 3) / 96.;
	p += pow(b, 5) * (d / 512. - 5 / 1536.);
	p -= pow(b, 7) * (d / 1536. - 29 / 23040.);

	return p;
}

double strong_coupling_pl_su3(double b, unsigned d){
	double p = b / 2.;
	p += pow(b, 2) / 8.;
	p -= pow(b, 4) * 5 / 384.;
	p += pow(b, 5) * (d / 1296. - 113 / 20736.);
	p += pow(b, 6) * (d * 7 / 5184. - 931 / 414720.);
	p += pow(b, 7) * (d * 5 / 5184. - 1069 / 829440.);

	return p;
}
