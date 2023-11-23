! gfortran test_matmul.f90 helpers.f90 -o test_matmul.x -fopenmp -llapack -lopenblas
!
! test_matmul.f90  
!
program main
    use helpers
    use omp_lib
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none
    real(dp), allocatable, dimension(:,:) :: A, B, C
    real(dp), allocatable, dimension(:) :: x, y

    complex(dp), allocatable, dimension(:,:) :: D, E, F
    complex(dp), allocatable, dimension(:) :: u, v

    integer :: ndim, i
    real(dp) :: res
    complex(dp) :: cres
    integer :: ci
    character (len=1024) :: arg
    
    ndim = 1000
    ci = 1
    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-n") then
            call get_command_argument(ci+1,arg)
            read(arg, '(i32)') ndim
            ci = ci + 2
        else
            exit
        end if
        
    end do

    allocate(A(ndim,ndim))
    allocate(B(ndim,ndim))
    allocate(C(ndim,ndim))
    allocate(x(ndim))
    allocate(y(ndim))  
    
    allocate(D(ndim,ndim))
    allocate(E(ndim,ndim))
    allocate(F(ndim,ndim))
    allocate(u(ndim))
    allocate(v(ndim)) 

    call random_number(A)
    call random_number(x)
    B = 0.0d0
    do i=1,ndim
        B(i,i) = 1.0d0
    end do

    D = cmplx(A, 0.0d0, dp)
    E = cmplx(B, 0.0d0, dp)
    u = cmplx(x, 0.0d0, dp)

    ! testing dmul_mm
    print*, "Testing: 'dmul_mm'...... "
    call dmul_mm(A, B, C)
    ! print*, A
    ! print*, B
    ! print*, C
    print*, "Passed test: dmul_mm"

    !testing dmul_mv
    print*, "Testing: 'dmul_mv'...... "
    call dmul_mv(B, x, y)
    ! print*, x
    ! print*, B
    ! print*, y
    print*, "Passed test: dmul_mv"
    
    !testing ddot
    print*, "Testing: 'ddot'...... "
    call dmul_ddot(x, y, res)
    print*, "Passed test: ddot"


    ! testing zmul_mm
    print*, "Testing: 'zmul_mm'...... "
    call zmul_mm(D,E,F)
    ! print*, D
    ! print*, E
    ! print*, F
    print*, "Passed test: zmul_mm"

    !testing zmul_mv
    print*, "Testing: 'zmul_mv'...... "
    call zmul_mv(E, u, v)
    ! print*, u
    ! print*, E
    ! print*, v
    print*, "Passed test: zmul_mv"
    
    !testing ddot
    print*, "Testing: 'zdotc'...... "
    call zmul_zdotc(u, v, cres)
    print*, "Passed test: zdotc"

end program main