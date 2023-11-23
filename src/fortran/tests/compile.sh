compiler='ifort'

intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi

$compiler test_matmul.f90 ../helpers.f90 -o test_matmul.x $flags 
$compiler test_ho_dvr.f90 ../helpers.f90 ../integrators.f90 -o test_ho_dvr.x $flags