#ifndef GAUGE_FFTW3
#define GAUGE_FFTW3
#include <fftw3.h>
#endif

#ifndef GAUGE_HARMONIC
#define GAUGE_HARMONIC

void fill_harm_mat(double complex *u, double beta, unsigned nn2, unsigned nl, unsigned i, gauge_flags *mode);

void sample_harmonic(double *x, double complex *xc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);
void sample_fourier_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);

double energy_fourier_momenta(double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);

void evolve_fields(double complex *xc, double beta, unsigned ns, unsigned nn2, double h, const fftw_plan *fft, gauge_flags *mode);

#endif
