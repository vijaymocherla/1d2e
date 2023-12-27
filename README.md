# 1d2e

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
- `imag_tprop.f90` contains the module `imag_tprop` with routines for imaginary-time propagation. `itp_on_the_fly` uses on-the-fly calculation of matrix elements and `itp_sparse` uses sparse matrix-vector multiplication routines.
- `real_tprop.f90` contains the module `real_tprop` with sparse matrix-vector multiplication routines for real-time propagation  (use `rtp_sparse`).
- `hartreefock.f90` contains the `hartreefock` module, with an implementations of self-consistent field (SCF) routine for an 2e- DVR Hartree-Fock calculation.
- `intracules.f90` contains a the `intracules` module, with routines that compute different intracules for a given wavefunction input.

## Tests

`tests/` contains:

- `test_1e_dvr.f90` is a test for 1e- single well calculation.
- `test_2e_dvr.f90` is a test for exact 2e- DVR with imaginary time propagation. 
- `test_ho_dvr.f90` is a test for subroutines in `integrators` using quantum harmonic oscillator problem.  
- `test_matmul.f90` is a test for the wrapper subroutines in `blas_wrappers.f90`.
- `test_hartreefock.f90` is a test for hartree-fock 2e- DVR calculation.
- `test_intracules.f90` is a test that computes intracules for given wavefunction.
- `test_sparse_blas.f90` is test for sparse BLAS routines used in the modules `imag_tprop` and `real_tprop`.
- `test_sparse_itp.f90` is test for imaginary-time propagation using sparse routine `itp_sparse` in module `imag_tprop`.
- `test_sparse_rtp.f90` is test for real-time propagation using sparse routine `rtp_sparse` in module `real_tprop`.