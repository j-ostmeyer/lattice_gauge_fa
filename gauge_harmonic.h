#ifndef GAUGE_FFTW3
#define GAUGE_FFTW3
#include <fftw3.h>
#endif

#ifndef GAUGE_HARMONIC
#define GAUGE_HARMONIC

void sample_id(double complex *u, unsigned ns, unsigned nn2, unsigned gd);

double fill_harm_mat(double complex *m, unsigned nn2, unsigned nl, unsigned i, gauge_flags *mode);

void mat_mul_basis(double complex *x, double complex *y, double complex *m, unsigned nn2, unsigned ng, unsigned j, int dagger);

void sample_fourier_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);
double energy_fourier_momenta(double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);
void get_fourier_x_dot(double complex *pc, double beta, unsigned ns, unsigned nn2, double h, const fftw_plan *fft, gauge_flags *mode);

#endif
