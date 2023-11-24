module helpers
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    
    private

    public zexp_array, dexp_array

    contains

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