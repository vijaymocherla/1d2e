program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use blas_wrappers
    implicit none

    real(dp), allocatable, dimension(:) :: A
    integer, allocatable, dimension(:) :: arow , acol
    real(dp), allocatable, dimension(:) :: X,Y
    integer :: n, i, j, k, nnz

    n = 3
    nnz = 3*(n-2)+4
    allocate(A(nnz), acol(nnz))
    allocate(arow(n+1))
    allocate(X(n), Y(n))

    k = 3
    A(1) = 2.0
    A(2)  = 1.0
    acol(1) = 1
    acol(2) = 2
    arow(1) = 1
    do i=2,n-1
        arow(i) = k
        A(k) = 1.0
        A(k+1) = 2.0
        A(k+2) = 1.0
        acol(k) = i-1
        acol(k+1) = i
        acol(k+2) = i+1
        k = k + 3
    end do
    arow(n) = k
    A(nnz-1) = 1.0
    A(nnz) = 2.0
    acol(nnz-1) = n-1
    acol(nnz) = n 
    arow(n+1) = nnz+1
    print*, k+1
    X = 1.0d0
    Y = 0.0d0
    print*, A
    print*, "arow"
    print*, arow
    print*, "acol"
    print*, acol
    print*, ""
    call csr_dmul_mv(A, arow, acol, X, Y)
    print*, Y


end program main