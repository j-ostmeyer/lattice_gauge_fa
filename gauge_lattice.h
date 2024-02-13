#ifndef GAUGE_LATTICE
#define GAUGE_LATTICE

unsigned *construct_lattice(FILE *input, unsigned *ns, unsigned *nn, gauge_flags *mode);
unsigned *construct_cubic(FILE *input, unsigned *ns, unsigned *nn, gauge_flags *mode);

#endif
