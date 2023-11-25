program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use omp_lib
    implicit none

    real(dp), allocatable, dimension(:,:) :: a
    real(dp), allocatable, dimension(:) :: x, y

    integer :: ndim
    integer :: i,j

    ndim = 500

    allocate(a(ndim, ndim))
    allocate(x(ndim), y(ndim))
    
    call random_number(a)
    call random_number(x)

    print*, "started loop"
    !$omp parallel 
    !$omp do    
        do i=1,ndim
            y(i) = 0.0d0
            do j=1,ndim
                y(i) = y(i) + a(i,j)*x(i)
                print*, i, j 
            end do
        end do
    !$omp end do
    !$omp end parallel
    print*, "ended loop"
    
    print*, "OMP! works"
end program main