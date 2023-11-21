flags="-fopenmp -llapack -lopenblas"

gfortran test_matmul.f90 helpers.f90 -o test_matmul.x $flags 
gfortran test_ho_dvr.f90 helpers.f90 integrators.f90 -o test_ho_dvr.x $flags