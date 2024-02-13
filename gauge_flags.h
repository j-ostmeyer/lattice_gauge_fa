#ifndef GAUGE_FLAGS
#define GAUGE_FLAGS

typedef struct GaugeFlags{
	int gauge_group;
	unsigned gauge_dim;
	unsigned num_gen;
	unsigned length_cube;
	unsigned space_dim;
	int no_fourier_acc;
	double complex *zdummy;
	double *ddummy;
	int *idummy;
} gauge_flags;

void set_flags(FILE *in, gauge_flags *mode);

#endif
