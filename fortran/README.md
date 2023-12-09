# Wigner

Fortran code for running 1-dimensional 2e- calculations.

> Note: `CMAKE_FC_COMPILER` is set to `ifort`.

To compile create a `build` directory and use CMake to compile scripts.
```sh
$ mkdir build & cd ./build;
$ cmake ..
$ make
```

Run calculations/tests using executables in `build/tests/`. 


## Source Files

`src/` contains:

- `blas_wrappers.f90` contains subroutines that wrap BLAS functions for vector-vector, matrix-vector and matrix-matrix operations (also contains Sparse BLAS wrappers).
- `lapack_wrappers.f90` contains subroutines that wrap LAPACK functions.
- `helpers.f90` contains  array operation routines parallelized with OpenMP and other miscellaneous helper subroutines.
- `two_electron_dvr.f90` contains the `two_electron_dvr` module with functions to evaluate matrix elements of 1e- and 2e- DVR hamiltonians.
- `rungekutta.f90` contains the `rungekutta` module, with simple implementations of real-time and imaginary-time propagation (these routines may not scale well for large systems).
- `itprop.f90` contains the module `itp` with routines for imaginary-time propagation. `itp_on_the_fly` uses on-the-fly calculation of matrix elements and `itp_sparse` uses sparse matrix-vector multiplication routines.

## Tests

`tests/` contains:

- `test_1e_dvr.f90` is a test for 1e- single well calculation.
- `test_2e_dvr.f90` is a test for 2e- DVR with imaginary time propagation. 
- `test_ho_dvr.f90` is a test for subroutines in `integrators` using quantum harmonic oscillator problem.  
- `test_matmul.f90` is a test for the wrapper subroutines in `blas_wrappers.f90`.