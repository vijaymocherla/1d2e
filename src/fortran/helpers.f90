module helpers
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    private

    public zmul_mv, zmul_mm, zmul_zdotc, &
           dmul_mv, dmul_mm, dmul_ddot, &
           zexp_array, dexp_array, & 
           eigsh 

    contains

    ! wrapper for ZGEMV
    subroutine zmul_mv(A, X, Y)
        complex(dp), dimension(:,:), intent(in) :: A
        complex(dp), dimension(:), intent(in) :: X
        complex(dp), dimension(:), intent(out) :: Y
        integer :: m, n, incx, incy, lda
        character :: transa
        complex(dp) :: alpha, beta
        lda = size(A,1)
        incx = 1
        incy = 1
        m = size(A,dim=1)
        n = size(A,dim=2)
        alpha = (1.0d0, 0.0d0)
        beta = (1.0d0, 0.0d0)
        transa = 'N'
        call zgemv(transa, m, n, alpha, A, lda, X, incx, alpha, Y, incy)
    end subroutine  

    ! wrapper for ZGEMM
    subroutine zmul_mm(A, B, C)
        complex(dp), dimension(:, :), intent(in) :: A, B
        complex(dp), dimension(:, :), intent(out) :: C
        integer :: l, m, n, lda, ldb, ldc
        character :: transa, transb, transc
        complex(dp) :: alpha, beta 
        alpha=(1.0d0, 0.0d0)
        beta=(1.0d0, 0.0d0)
        transa = 'N'
        transb = 'N'
        transc = 'N'
        lda = size(A,1)
        ldb = size(B,1)
        ldc = size(C,1)
        l = size(A,1)
        n = size(B,1)
        m = size(C,1)
        call zgemm(transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine
    
    ! wrapper for ZDOTU
    subroutine zmul_zdotu(X, Y, res)
        complex(dp), dimension(:), intent(in) :: X, Y 
        complex(dp), intent(out) :: res
        complex(dp), external :: zdotu
        integer :: incx, incy, ndim
        ndim = size(X, dim=1)
        incx = 1
        incy = 1
        res = zdotu(ndim, x, incx, y, incy)
    end subroutine  

    ! wrapper for ZDOTC
    subroutine zmul_zdotc(X, Y, res)
        complex(dp), dimension(:), intent(in) :: X, Y 
        complex(dp), intent(out) :: res
        complex(dp), external :: zdotc
        integer :: incx, incy, ndim
        ndim = size(X)
        incx = 1
        incy = 1
        res = zdotc(ndim, x, incx, y, incy)
    end subroutine      

    ! wrapper for DGEMV
    subroutine dmul_mv(A, X, Y)
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(:), intent(in) :: X
        real(dp), dimension(:), intent(out) :: Y
        integer :: m, n, incx, incy, lda, ndim
        real(dp) :: alpha, beta
        character :: transa
        ndim = size(A, dim=1)
        lda = size(A,dim=1)        
        incx = 1
        incy = 1
        m = size(A,dim=1)
        n = size(x,dim=1)
        alpha = 1.0d0
        beta = 1.0d0
        transa = 'N'
        call dgemv(transa, m, n, alpha, A, lda, X, incx, beta, Y, incy)
    end subroutine  

    ! wrapper for DGEMM
    subroutine dmul_mm(A, B, C)
        real(dp), dimension(:, :), intent(in) :: A, B
        real(dp), dimension(:, :), intent(out) :: C
        integer :: l, m, n, lda, ldb, ldc, ndim
        character :: transa, transb, transc
        real(dp) :: alpha, beta 
        alpha = 1.0d0
        beta = 1.0d0
        transa = 'N'
        transb = 'N'
        transc = 'N'
        ndim = size(A,dim=1) 
        lda = size(A,dim=1)
        ldb = size(B,dim=1)
        ldc = size(C,dim=1)
        l = ndim
        n = ndim
        m = ndim
        call dgemm(transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

    ! wrapper for DDOT
    subroutine dmul_ddot(X, Y, res)
        real(dp), dimension(:), intent(in) :: X, Y 
        real(dp), intent(out) :: res
        real(dp), external :: ddot
        integer :: incx, incy, ndim
        ndim = size(X, dim=1)
        incx = 1
        incy = 1
        res = ddot(ndim, x, incx, y, incy)
    end subroutine

    ! wrapper for DSYEV
    subroutine eigsh(A, vals, vecs)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(out) :: vals
        real(dp), dimension(:, :), intent(out) :: vecs
        integer :: info, ndim
        real(dp), allocatable :: work(:)
        character, parameter :: jobz="V", uplo="U"
        ndim = size(A, dim=1) 
        allocate(work(3*ndim))
        vecs = A
        call dsyev(jobz, uplo, ndim, vecs, ndim, vals, work, 3*ndim, info)
    end subroutine eigsh

    ! Returns expm(1j*alpha*x)
    subroutine zexp_array(alpha, x, exp_x)
        real(dp), intent(in) :: alpha
        complex(dp), dimension(:), intent(in) :: x
        complex(dp), dimension(:,:), intent(out) :: exp_x
        integer :: i, ndim
        complex(dp) :: tmp
        ndim = size(x)
        do i = 1, ndim
            tmp = dcmplx(0.0d0, alpha)*x(i)
            exp_x(i,i) = exp(tmp) 
        end do
    end subroutine

    ! Returns expm(alpha*x)
    subroutine dexp_array(alpha, x, exp_x)
        real(dp), intent(in) :: alpha
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:,:), intent(out) :: exp_x
        integer :: i, ndim
        real(dp) :: tmp
        ndim = size(x)
        do i = 1, ndim
            tmp =  alpha*x(i)
            exp_x(i,i) = exp(tmp) 
        end do
    end subroutine


end module helpers