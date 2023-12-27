module intracules
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use omp_lib
    implicit none
    
    private

    public :: comp_wigner_intracule

    real(dp), parameter, public :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 

    contains

    subroutine comp_wigner_intracule(psi, x0, wigner_intracule)
        implicit none
        complex(dp), allocatable, intent(in) :: psi(:) 
        real(dp), allocatable, intent(out) :: wigner_intracule (:)
        real(dp), intent(in) :: x0 
        integer :: ndim, n
        integer :: u, v
        integer :: i, j, m
        integer :: a, b
        real(dp) :: dx, L
        complex(dp) :: t1,t2
        ndim = size(psi)
        n = int(sqrt(real(ndim,kind=dp)))
        L = x0/4
        dx = 2*x0/(n-1)
        m = int(2*L/dx) + 1
        allocate(wigner_intracule(m*m))
        !$omp parallel
        !$omp do
        do u=1,m
            do v=1,m
                t1 = 0.0d0
                do i=1,m
                    t2 = 0.0d0
                    do j=1,m
                        a = i+j-m+(n/2)
                        b = j-i+(n/2)
                        t2 = t2 + dconjg(psi(a*n+b+u+1))*psi(b*n+a+u+1)
                        t2 = t2 + dconjg(psi(a*n+b-u+1))*psi(b*n+a-u+1)
                    end do 
                    t2 = t2*cos(2*v*dx*(i*dx-L))*2/(pi_dp)
                    t1 = t1 + t2
                end do
                wigner_intracule(m*(u-1)+v) = real(t1,kind=8)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine comp_wigner_intracule


end module intracules