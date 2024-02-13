#ifndef GAUGE_PLAQUETTE
#define GAUGE_PLAQUETTE

double plaquette_av(double *x, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double plaquette_fields(double *x, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);

void staple_fields(double *x, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);
void sum_of_plaquettes(double *x, double complex *pl, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, int direction, gauge_flags *mode);
void sum_of_staples(double *x, double complex *st, double complex *u, double complex *g, double *ev, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, gauge_flags *mode);

double strong_coupling_plaquette(double beta, gauge_flags *mode);
double strong_coupling_pl_u1(double b, unsigned d);
double strong_coupling_pl_su2(double b, unsigned d);
double strong_coupling_pl_su3(double b, unsigned d);

#endif
