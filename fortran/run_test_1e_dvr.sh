compiler='ifort'

intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi

cd tests/

$compiler -c test_1e_dvr.f90 ../src/exact_2e_dvr.f90 ../src/blas_wrappers.f90 ../src/lapack_wrappers.f90 $flags -heap-arrays
$compiler test_1e_dvr.o lapack_wrappers.o blas_wrappers.o exact_2e_dvr.o -o test_1e_dvr.x $flags -heap-arrays
./test_1e_dvr.x 


