#ifndef GAUGE_PLAQUETTE
#define GAUGE_PLAQUETTE

double plaquette_av(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double complex plaquette_tr(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);
void plaquette_mat(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, int direction, gauge_flags *mode);

void sum_of_plaquettes(double complex *u, double complex *pl, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, int direction, gauge_flags *mode);
void sum_of_staples(double complex *u, double complex *st, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, gauge_flags *mode);
void staple_fields(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);
void copy_staple(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);

double topo_charge(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double topo_charge_2d(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double topo_charge_4d(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode);
double complex clover_field_tr(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);
void clover_mat(double complex *u, double complex *cl, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode);

double wilson_loop_av(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned a, unsigned b, gauge_flags *mode);
double complex wilson_loop_tr(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, unsigned a, unsigned b, gauge_flags *mode);
unsigned prod_path(double complex *u, double complex *path, double complex *prod, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned a, gauge_flags *mode);

double strong_coupling_plaquette(double beta, gauge_flags *mode);
double strong_coupling_pl_u1(double b, unsigned d);
double strong_coupling_pl_pade_u1(double b, long d);
double strong_coupling_pl_su2(double b, unsigned d);
double strong_coupling_pl_pade_su2(double b, long d);
double strong_coupling_pl_su3(double b, unsigned d);
double strong_coupling_pl_pade_su3(double b, double d);

#endif
