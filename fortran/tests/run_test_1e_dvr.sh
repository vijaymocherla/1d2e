compiler='ifort'

intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

modules='../src/'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi

cd tests/

$compiler -c test_1e_dvr.f90 $flags -heap-arrays
$compiler test_1e_dvr.o ../src/* $flags -heap-arrays
# ./test_1e_dvr.x 


