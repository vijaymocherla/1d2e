module blas_wrappers
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    private

    public zmul_mv, zmul_mm, zmul_zdotc, &
           dmul_mv, dmul_mm, dmul_ddot
    contains

    ! wrapper for ZGEMV
    subroutine zmul_mv(A, X, Y, alpha)
        complex(dp), intent(in) :: alpha
        complex(dp), dimension(:,:), intent(in) :: A
        complex(dp), dimension(:), intent(in)   :: X
        complex(dp), dimension(:), intent(out)  :: Y
        integer :: m, n, incx, incy, lda
        character :: transa
        complex(dp) ::  beta
        lda = size(A,1)
        incx = 1
        incy = 1
        m = size(A,dim=1)
        n = size(A,dim=2)
        beta = (0.0d0, 0.0d0)
        transa = 'N'
        call zgemv(transa, m, n, alpha, A, lda, X, incx, beta, Y, incy)
    end subroutine  

    ! wrapper for ZGEMM
    subroutine zmul_mm(A, B, C)
        complex(dp), dimension(:,:), intent(in) :: A, B
        complex(dp), dimension(:,:), intent(out) :: C
        integer :: l, m, n, lda, ldb, ldc
        character :: transa, transb, transc
        complex(dp) :: alpha, beta 
        alpha=(1.0d0, 0.0d0)
        beta=(0.0d0, 0.0d0)
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
    subroutine dmul_mv(A, X, Y, alpha)
        real(dp), intent(in) :: alpha
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(:),   intent(in) :: X
        real(dp), dimension(:),  intent(out) :: Y
        integer :: m, n, incx, incy, lda, ndim
        real(dp) :: beta
        character :: transa
        ndim = size(A, dim=1)
        lda = size(A,dim=1)        
        incx = 1
        incy = 1
        m = size(A,dim=1)
        n = size(x,dim=1)
        beta = 0.0d0
        transa = 'N'
        call dgemv(transa, m, n, alpha, A, lda, X, incx, beta, Y, incy)
    end subroutine  

    ! wrapper for DGEMM
    subroutine dmul_mm(A, B, C)
        real(dp), dimension(:,:), intent(in)  :: A, B
        real(dp), dimension(:,:), intent(out) :: C
        integer :: l, m, n, lda, ldb, ldc, ndim
        character :: transa, transb, transc
        real(dp) :: alpha, beta 
        alpha = 1.0d0
        beta = 0.0d0
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
        real(dp), allocatable, dimension(:), intent(in) :: X, Y 
        real(dp), intent(out) :: res
        real(dp), external :: ddot
        integer :: incx, incy, ndim
        ndim = size(X, dim=1)
        incx = 1
        incy = 1
        res = ddot(ndim, x, incx, y, incy)
    end subroutine


end module blas_wrappers