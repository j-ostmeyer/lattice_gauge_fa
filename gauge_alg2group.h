#ifndef GAUGE_ALG2GROUP
#define GAUGE_ALG2GROUP

#define ISQ3 0.5773502691896258

void alg2group(double *x, double complex *u, double complex *g, double *ev, int dagger, gauge_flags *mode);
void get_alg_u1(double *x, double complex *g);
void get_alg_su2(double *x, double complex *g);
void get_alg_su3(double *x, double complex *g);

void project_tr_lambda(double complex *u, double *x, gauge_flags *mode);
void project_tr_u1(double complex *u, double *x);
void project_tr_su2(double complex *u, double *x);
void project_tr_su3(double complex *u, double *x);

#endif
