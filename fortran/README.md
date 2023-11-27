# Wigner

Fortran code for running 1-dimensional 2e- calculations.

> Note: If using intel compilers set the option `compiler=ifort` in the following bash scripts .

- Use `run_test_1e_dvr.sh` to compile and run Harmonic oscillator test for `two_electron_dvr` module.
- Use `run_test_2e_dvr.sh` to compile and run 2e- calculation test for `two_electron_dvr` module.


## Source Files

`src/` contains:

- `blas_wrappers.f90` contains subroutines that wrap BLAS functions for vector-vector, matrix-vector and matrix-matrix operations.
- `lapack_wrappers.f90` contains subroutines that wrap LAPACK functions.
- `helpers.f90` contains other miscellaneous helper subroutines.
- `exact_2e_dvr.f90` contains the `two_electron_dvr` module, also contains an imaginary-time propagation code paralleised with `OpenMP`.
- `integrators.f90` contains the `rungekutta` module, with implementations of real-time and imaginary-time propagation.

## Tests

`tests/` contains:

- `test_1e_dvr.f90` is a test for 1e- calculation using the module `two_electron_dvr`.
- `test_2e_dvr.f90` is a test for 2e- calculation using the module `two_electron_dvr`.
- `test_ho_dvr.f90` is a test for subroutines in `integrators` using quantum harmonic oscillator problem.  
- `test_matmul.f90` is a test for the wrapper subroutines in `blas_wrappers.f90`.

