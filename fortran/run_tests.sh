compiler='ifort'

intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi

cd tests/

# $compiler -c test_matmul.f90 -I../bin
# $compiler -c test_ho_dvr.f90 -I../bin

$compiler  test_matmul.f90 ../src/* -o test_matmul.x $flags 
./test_matmul.x

$compiler  test_ho_dvr.f90 ../src/* -o test_ho_dvr.x $flags -heap-arrays
./test_ho_dvr.x 


