include_directories(${CMAKE_BINARY_DIR}/tests/)
include_directories(${MKL_INCLUDE_DIR})

file(GLOB_RECURSE tests ${CMAKE_BINARY_DIR}/../src/*.f90)

set(CMAKE_C_COMPILER icx)
set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_FC_COMPILER ifort)

find_package(MKL REQUIRED)
find_package(OpenMP REQUIRED)
set(MKL_OMP_LIB_PATHS "")
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARIES})
# find_package(LAPACK REQUIRED)
# find_package(BLAS REQUIRED)


set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_FC_FLAGS "-O3 -qmkl")


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



target_link_libraries(test_ho_dvr ${MKL_LIBRARIES})
target_link_libraries(test_1e_dvr ${MKL_LIBRARIES})
target_link_libraries(test_2e_dvr ${MKL_LIBRARIES})

