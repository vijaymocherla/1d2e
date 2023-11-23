compiler='ifort'

intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi

$compiler helpers.f90 -c $flags
$compiler integrators.f90 -c $flags