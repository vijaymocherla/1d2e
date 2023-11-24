module integrators
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    private

    public rungekutta, imagtp
    
    contains

    subroutine rungekutta(func, y0, ti, tf, dt, print_nstep, outfile)
        use blas_wrappers
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        complex(dp), intent(in) :: func(:,:)
        complex(dp), intent(in) :: y0(:)
        real(dp), intent(in) :: ti, tf, dt
        integer, intent(in) :: print_nstep
        character (len=*), intent(in) :: outfile
       
        complex(dp), allocatable :: yi(:)
        complex(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        complex(dp), allocatable, dimension(:) :: kt

        
        complex(dp) :: expt, norm, autocorr
        real(dp) :: t
        integer :: n, m
        integer :: tstep
        
        n = size(func, dim=2)
        m = size(y0)
        ! some consistency checks
        if (n .ne. m) then
            print*, "!!ERROR: Dimensions of func and y0 don't match."
            stop
        end if

        allocate(k1(n), k2(n), k3(n), k4(n))
        allocate(kt(n))
        allocate(yi(n))
    
        t = ti
        tstep = 0
        yi = y0
        open(100, file=trim(outfile))
        ! writing for t0
        call zmul_zdotc(yi, yi, norm)
        ! expectation values
        call zmul_mv(func, yi, kt)
        call zmul_zdotc(yi,kt, expt)
        ! autocorr
        call zmul_zdotc(y0, yi, autocorr)
        write(100,'(4f32.16)') t, abs(norm), abs(autocorr), abs(expt)  

        do while (t < tf)
            call zmul_mv(func, yi, k1)
            call zmul_mv(func, (yi + (dt*0.5d0)*k1), k2)
            call zmul_mv(func, (yi + (dt*0.5d0)*k2), k3)
            call zmul_mv(func, (yi + (dt*1.0d0)*k3), k4)
            yi = yi + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
            t = t + dt
            tstep = tstep + 1
            if (modulo(tstep, print_nstep) == 0) then
                ! norm
                call zmul_zdotc(yi, yi, norm)
                ! intermediate normalisation only for itp
                yi = 1/(norm)**(0.5d0) * yi
                call zmul_zdotc(yi, yi, norm)
                ! expectation values
                call zmul_mv(func, yi, kt)
                call zmul_zdotc(yi,kt, expt)
                ! autocorr
                call zmul_zdotc(y0, yi, autocorr)
                write(100,'(4f32.16)') t, abs(norm), abs(autocorr), abs(expt)  
            end if
        end do
        close(100)
    end subroutine rungekutta

    ! subroutine for imaginary time propagation
    subroutine imagtp(H, psi, dt, nstep)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        use blas_wrappers
        implicit none
        real(dp), allocatable, intent(in) :: H(:,:)
        real(dp), allocatable, intent(inout) :: psi(:)
        integer, intent(in) :: nstep
        real(dp), intent(inout) :: dt
        complex(dp), allocatable :: func(:,:)
        complex(dp), allocatable :: y0(:), y(:)
        complex(dp) :: expt
        real(dp)  :: ti  ! initial time
        real(dp)  :: tf  ! final time
        integer :: n
    
        n = size(H, dim=1)
        allocate(y0(n), y(n))
        allocate(func(n,n))
        ! converting fs to au
        ti = 0.0d0
        tf = nstep*dt

        ! real time prop
        ! func = cmplx(0.00d0, -H, dp)
        ! imaginary time prop 
        func = cmplx(-H, 0.0d0)
        y0 = cmplx(psi, 0.0d0, dp)

        
        call rungekutta(func, y0, ti, tf, dt, print_nstep=1, outfile='itp.out')
    end subroutine

end module integrators