! gfortran test_ho_dvr.f90 helpers.f90 integrators.f90 -o test_ho_dvr.x -llapack -lopenblas -fopenmp 
!
! test_ho_dvr.f90 
!
program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use blas_wrappers, only: dmul_ddot, dmul_mv
    use lapack_wrappers, only: eigsh
    use integrators, only: rungekutta, imagtp
    implicit none
    real(dp), allocatable, dimension(:,:) :: H
    real(dp), allocatable, dimension(:,:) :: T
    real(dp), allocatable, dimension(:,:) :: V
    real(dp), allocatable, dimension(:)   :: x
    real(dp), allocatable, dimension(:)   :: psi0
    real(dp), allocatable, dimension(:)   :: evals
    real(dp), allocatable, dimension(:,:) :: evecs
    
    real(dp), parameter  :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 

    real(dp)  :: x0  ! grid size (-x0, x0)
    integer   :: n   ! number of grid points 
    real(dp)  :: dx  ! grid spacing 
    real(dp)  :: m   ! mass of the electron
    real(dp)  :: dt
    

    integer   :: i,j
    real(dp)  :: norm
    real(dp)  :: expt


    x0 = 10.0d0
    n = 100
    dx = 2*x0/n
    m = 1.0d0
    allocate(T(n,n), V(n,n), H(n,n))
    allocate(x(n), psi0(n))
    H = 0.0d0
    T = 0.0d0
    V = 0.0d0
    
    do i=1,n
        x(i) = -x0 + (i-1)*dx
    end do

    ! Kinetic operator
    do i = 1, n
        T(i,i) = (-1)**(i-i) * (1.0d0 / (2.0d0 * m * dx**2)) * pi_dp**2 / 3.0d0
        do j = i+1, n
            T(i,j) = (-1)**(i-j) * (1.0d0 / (2.0d0 * m * dx**2)) * 2.0d0 / (i-j)**2
            T(j,i) = T(i,j)
        enddo
    enddo

    ! Potential operator
    do i=1,n
        V(i,i) = 0.50d0*x(i)**2
    end do

    H = T + V
    deallocate(T,V)
    

    ! digonalising
    allocate(evals(n), evecs(n,n))
    call eigsh(H, evals, evecs)


    open(100,file='evals.txt')
        do i=1,n
            write(100,*) i, evals(i)
        end do
    close(100)


    psi0 = 0.0d0
    dt  = 0.01d0
    call random_number(psi0)
    call dmul_ddot(psi0, psi0, norm)
    psi0 = 1/(norm)**(0.5d0) * psi0

    call imagtp(H, psi0, dt, 200)


end program main

