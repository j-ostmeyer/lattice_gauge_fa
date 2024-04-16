# Lattice Gauge Theory with HMC and Fourier Acceleration

by Johann Ostmeyer

An implementation of pure lattice gauge theory simulations with Hybrid Monte Carlo (HMC) and the optional feature of Fourier acceleration.

This repository contains the code required to reproduce the results presented in the section on lattice gauge theory of *"Minimal Autocorrelation in Hybrid Monte Carlo simulations with Exact Fourier Acceleration"*, [arXiv:2404.09723 [hep-lat]](https://arxiv.org/abs/2404.09723).

For questions concerning the code or data contact [ostmeyer@hiskp.uni-bonn.de](mailto:ostmeyer@hiskp.uni-bonn.de).

## Implementation

The entire simulations are implemented in `C`.

### Compile
The code depends on `Lapacke` and `FFTW3`. Make sure you have recent installations available before you proceed.

After adjusting the Makefile according to your needs (possibly switching to another compiler), simply type `make` to compile all the `C` files.

### Run Code
Adjust the `input.txt` file according to the gauge group, dimensionality and method you want to simulate. Alternatively, create a new file with the relevant input parameters. Then execute:
```
./gauge_main input.txt results.csv
```
As the name suggests, all results will be written to `results.csv` (or any other file provided as a second argument).

### Results
The results come in at least 7 columns:

1. Average plaquette
2. deviation of plaquette from strong coupling expansion
3. topological charge (only tested for U(1) in 2D)
4. squared topological charge
5. HMC energy
6. Acceptance (acc=1, rej=0)
7. Boltzmann weight `exp(-(E-E_0))`\
... any number of Wilson loops (passed sanity checks, but not tested extensively) from $2\times2$ to $l\times l$, $l$ specified on input

## Data

Most of the simulations in paper can be reproduced very quickly, nonetheless the resulting data will gladly be provided upon request.

## Stable Releases

`v1.1.1` initial publication with the preprint, compatible with arXiv version v1.
