#ifndef GAUGE_FFTW3
#define GAUGE_FFTW3
#include <fftw3.h>
#endif

#ifndef GAUGE_HMC
#define GAUGE_HMC

void gauge_force(double complex *u, double *p_dot, double beta, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);

double energy_momenta(double *p, double complex *pc, double beta, unsigned ns, unsigned nn2, const fftw_plan *fft, gauge_flags *mode);
double hamilton(double energyP, double pl_av, double beta, unsigned ns, gauge_flags *mode);

void update_u(double complex *u, double *p, double complex *pc, double *x_dot, double beta, unsigned ns, unsigned nn, double h, const fftw_plan *fft, gauge_flags *mode);

void leap_frog(double complex *u, double *p, double complex *pc, double *p_dot, double *x_dot, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, const fftw_plan *fft, gauge_flags *mode);
short trajectory(double *p, double complex *pc, double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, int iteration, unsigned meas_freq, double *plaquettes_old, const fftw_plan *fft, FILE *res_out, gauge_flags *mode);

void run_hmc(double beta, unsigned *nnt, unsigned ns, unsigned nn, unsigned steps, double traj_length, unsigned therm, unsigned meas, unsigned meas_freq, const char *res_name, gauge_flags *mode);

#endif
