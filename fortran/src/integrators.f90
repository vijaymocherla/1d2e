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
    subroutine imagtp(H, psi0, dt, print_nstep, etol)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        use blas_wrappers
        implicit none
        real(dp), dimension(:,:), intent(in) :: H
        real(dp), dimension(:), intent(inout) :: psi0
        integer, intent(in) :: print_nstep
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: etol
        real(dp)  :: t  ! time
        integer :: n
        integer :: tstep
        real(dp), allocatable :: psi_i(:)
        real(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        real(dp), allocatable, dimension(:) :: kt
        real(dp) :: norm, autocorr, E0, Ei, dE
    
        n = size(H, dim=1)
        allocate(k1(n), k2(n), k3(n), k4(n))
        allocate(kt(n))
        allocate(psi_i(n))
        ! converting fs to au
        t = 0.0d0
        tstep = 0
        psi_i = psi0
        E0 = 0.0d0
        Ei = 0.0d0
        call dmul_mv(H, psi_i, kt)
        call dmul_ddot(kt, psi_i, E0)
        open(100, file='itp.out')
        do 
            call dmul_mv(-H, psi_i, k1)
            call dmul_mv(-H, (psi_i + (dt*0.5d0)*k1), k2)
            call dmul_mv(-H, (psi_i + (dt*0.5d0)*k2), k3)
            call dmul_mv(-H, (psi_i + (dt*1.0d0)*k3), k4)
            psi_i = psi_i + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
            t = t + dt
            tstep = tstep + 1
            if (modulo(tstep, print_nstep) == 0) then
                ! norm
                call dmul_ddot(psi_i, psi_i, norm)
                ! intermediate normalisation only for itp
                psi_i = 1/(norm)**(0.5d0) * psi_i
                call dmul_ddot(psi_i, psi_i, norm)
                ! autocorr
                call dmul_ddot(psi0, psi_i, autocorr)
                ! calculating energy
                call dmul_mv(H, psi_i, kt)
                call dmul_ddot(kt, psi_i, Ei)
                dE = abs(E0 - Ei)
                E0 = Ei
                write(100,'(4f32.16)') t, norm, autocorr, Ei 
                if (dE < etol) then
                    exit
                end if
            end if
        end do
        close(100)

    end subroutine

end module integrators