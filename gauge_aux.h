#ifndef GAUGE_AUX
#define GAUGE_AUX

#define NUM_RES 7

#ifndef M_PI
#define M_PI 3.141592653589793
#endif
#define M_2PI 6.283185307179586
#define HALF_M_PI 1.570796326794897

#define max(a,b)              \
	({                        \
	 __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b;       \
	 })

#define min(a,b)              \
	({                        \
	 __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b;       \
	 })

double norm(double complex x);
double vec_norm(double complex *x, unsigned n);

void print_vec(double *vec, unsigned n);
void print_zvec(double complex *vec, unsigned n);
void print_mat(double *m, unsigned long n);
void print_zmat(double complex *m, unsigned long n);
void fprint_results(FILE *out, double *res, unsigned n, unsigned r);

double average(double *x, unsigned n);
double average_sq(double *x, unsigned n);
double complex scalar_dot(double complex *x, double complex *y, unsigned d);

double reTr(double complex *u, unsigned d);
void construct_id(double complex *m, unsigned n);

void check_unitarity(double complex *u, unsigned ns, unsigned nn, gauge_flags *mode);

void copy_mat(double complex *x, double complex *y, unsigned n, int dagger);

void dagger(double complex *m, unsigned n);
void dagger_sym(double complex *m, unsigned n);
void dagger_asym(double complex *m, unsigned n);

void mat_mul(double complex *x, double complex *y, double complex *res, unsigned n);
void mat_prod(double complex *u, unsigned n, unsigned d);
double trace_prod(double complex *u, unsigned n, unsigned d);

#endif
