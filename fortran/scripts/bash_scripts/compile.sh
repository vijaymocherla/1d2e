compiler='ifort'
source='../src/'
intel_flags='-qmkl -qopenmp'
gnu_flags='-llapack -lopenblas -fopenmp'

if [ $compiler == 'ifort' ]; then
    flags=$intel_flags
elif [ $compiler == 'gfortran' ]; then
    flags=$gnu_flags
fi
cd bin/

$compiler $source'blas_wrappers.f90' -c $flags
$compiler $source'lapack_wrappers.f90' -c $flags
$compiler $source'helpers.f90' -c $flags
$compiler $source'itprop.f90' -c $flags

$compiler $source'two_electron_dvr.f90' -c $flags
$compiler $source'rungekutta.f90' -c $flags
