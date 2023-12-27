include_directories(${CMAKE_BINARY_DIR}/tests/)
include_directories(${MKL_INCLUDE_DIR})

file(GLOB_RECURSE tests ${CMAKE_BINARY_DIR}/../src/*.f90)
set(CMAKE_C_COMPILER icx)
set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_Fortran_COMPILER ifort)

find_package(MKL REQUIRED)
find_package(OpenMP REQUIRED)
set(MKL_OMP_LIB_PATHS "")
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARIES})
# find_package(LAPACK REQUIRED)
# find_package(BLAS REQUIRED)


set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_Fortran_FLAGS "-O3 -qmkl -qopenmp -funroll-loops -heap-arrays")

set(CMAKE_VERBOSE_MAKEFILE on)

add_executable(test_ho_dvr
            ${tests}
            test_ho_dvr.f90
) 

add_executable(test_1e_dvr
            ${tests}
            test_1e_dvr.f90
)               

add_executable(test_2e_dvr
            ${tests}
            test_2e_dvr.f90
)   

add_executable(test_sparse_itp
            ${tests}
            test_sparse_itp.f90
)   

add_executable(test_sparse_blas
            ${tests}
            test_sparse_blas.f90
) 

add_executable(test_sparse_rtp
            ${tests}
            test_sparse_rtp.f90
)               

add_executable(test_intracules
            ${tests}
            test_intracules.f90
)               

add_executable(test_hartreefock
            ${tests}
            test_hartreefock.f90
)   