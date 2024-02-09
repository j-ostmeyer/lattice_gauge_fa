CC = gcc
MYFILES=lattice flags aux random hmc
RNG=mt19937-64

CC += -std=gnu11 -march=native
MYFLAGS=-Wall -pedantic -O3
MYFLAGS += -Wno-unknown-pragmas
CC += -Wno-unused-result
MYFLAGS += -Wno-comment

MYLIBS=-lfftw3 -lm


findany = $(strip $(foreach W,$1,$(findstring $W,$2)))

ifneq (,$(call findany,login viz,$(HOSTNAME)))
#  #We are on DIaL
  CC = icc -std=c11
  #MYFLAGS += -qopenmp
  MYFLAGS += -no-multibyte-chars -mkl
  #MYFLAGS += -fopenmp
  MYFLAGS += -xCORE-AVX512 #-xHost #-ipo -fp-model fast=2
  MYFLAGS += -I$$FFTW3_DOUBLEINCLUDE
  MYFLAGS += -DMKL_LAPACKE
  MYLIBS += -lmkl_core -liomp5 -lmkl_intel_thread -lmkl_rt -lmkl_intel_lp64 #-lpthread
  MYLIBS += -L$$FFTW3_DOUBLELIB
  #MYFLAGS += -I$$LAPACKINCLUDE
  #MYLIBS += -L$$LAPACKLIB
  #MYLIBS += -L$$OPENBLASLIB -lopenblas
else
  #MYFLAGS += -g -Og -fsanitize=address -static-libasan
  MYLIBS += -llapacke -llapack -lblas
endif

all: gauge_main

# Take file identifiers and change empty prefix by organic_ % empty suffix by .o
# Compile all dependencies $^ to target $@
gauge_main: $(MYFILES:%=gauge_%.o) $(RNG:%=%.o) gauge_main.o
	$(CC) $(MYFLAGS) -o $@ $^ $(MYLIBS)

# For every object .o check for .c and .h file with same name
# Compile first dependency $< (the .c file)
%.o: %.c %.h Makefile
	$(CC) $(MYFLAGS) -c $<

print-%  : ; @echo $* = $($*)

clean:
	rm -f *.o *.so gauge_main

distclean:
	rm -f *.o *.so *.exe *~ tags gauge_main
