module helpers
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    use omp_lib
    use blas_wrappers, only:dmul_ddot
    implicit none
    
    private

    public zexp_array, dexp_array, &
           omp_normalize, omp_daxpy, omp_dotprod

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

    ! parallelised dot product of two vectors
    subroutine omp_dotprod(psi_i, psi_j, result)
        implicit none
        real(dp), dimension(:), intent(in)  :: psi_i, psi_j
        real(dp), intent(out) :: result
        integer :: u
        integer :: ndim
        result = 0.0d0
        ndim = size(psi_i)
        !$omp parallel private(u) shared(ndim, psi_i, psi_j, result)
        !$omp do reduction ( + : result )
        do u = 1, ndim
            result = result + psi_i(u) * psi_j(u)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine

    ! parallelised constant times a vector plus a vector ( with double precision) 
    subroutine omp_daxpy(a, x_array, y_array, z_array)
        implicit none
        real(dp), intent(in) :: a
        real(dp), dimension(:), intent(in)  :: x_array, y_array
        real(dp), dimension(:), intent(out) :: z_array
        integer :: i, ndim
        ndim = size(x_array)
        z_array = 0.0d0
        !$omp parallel private(i) shared(ndim, a, x_array, y_array, z_array)
        !$omp do
        do i=1,ndim
            z_array(i) = a*x_array(i) + y_array(i)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine

    ! parallelised vector normalisation
    subroutine omp_normalize(psi)
        real(dp), intent(inout) :: psi(:)
        real(dp) :: norm, inv_norm
        integer :: u, ndim
        ndim = size(psi)
        call omp_dotprod(psi, psi, norm)
        inv_norm = 1.0d0/(norm)**(0.5)
        !$omp parallel private(u) shared(ndim, inv_norm, psi)
        !$omp do
        do u=1,ndim
            psi(u) = inv_norm * psi(u)
        end do
        !$omp end do 
        !$omp end parallel
    end subroutine omp_normalize

    ! ! normalise an array
    ! subroutine normalize(psi)
    !     real(dp), intent(inout) :: psi(:)
    !     real(dp) :: norm, inv_norm
    !     integer :: ndim
    !     ndim = size(psi)
    !     call dmul_ddot(psi, psi, norm)
    !     inv_norm = 1.0d0/(norm)**(0.5)
    !     psi = inv_norm * psi
    ! end subroutine normalize    


end module helpers