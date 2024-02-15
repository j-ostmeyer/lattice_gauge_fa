#ifndef GAUGE_PLAQUETTE
#define GAUGE_PLAQUETTE

double plaquette_av(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double plaquette_fields(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);

void sum_of_plaquettes(double complex *u, double complex *pl, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, int direction, gauge_flags *mode);
void sum_of_staples(double complex *u, double complex *st, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, gauge_flags *mode);
void staple_fields(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);

double strong_coupling_plaquette(double beta, gauge_flags *mode);
double strong_coupling_pl_u1(double b, unsigned d);
double strong_coupling_pl_pade_u1(double b, double d);
double strong_coupling_pl_su2(double b, unsigned d);
double strong_coupling_pl_pade_su2(double b, double d);
double strong_coupling_pl_su3(double b, unsigned d);
double strong_coupling_pl_pade_su3(double b, double d);

#endif
