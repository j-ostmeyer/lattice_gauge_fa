#ifndef GAUGE_AUX
#define GAUGE_AUX

#define NUM_RES 5

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

void print_vec(double *vec, unsigned n);
void print_zvec(double complex *vec, unsigned n);
void print_mat(double *m, unsigned long n);
void print_zmat(double complex *m, unsigned long n);
void fprint_results(FILE *out, double *res, unsigned n, unsigned r);

double average(double *x, unsigned n);
double average_sq(double *x, unsigned n);
double scalar_dot(double *x, double *y, unsigned d);

double reTr(double complex *u, unsigned d);
void construct_id(double complex *m, unsigned n);

void dagger(double complex *m, unsigned n);
void dagger_sym(double complex *m, unsigned n);
void dagger_asym(double complex *m, unsigned n);

void mat_mul(double complex *x, double complex *y, double complex *res, unsigned n);
void mat_prod(double complex *u, unsigned n, unsigned d);
double trace_prod(double complex *u, unsigned n, unsigned d);

#endif
