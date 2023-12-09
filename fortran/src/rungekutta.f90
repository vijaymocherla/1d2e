module rungekutta
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    private

    public real_tprop, imag_tprop
    
    contains

    subroutine real_tprop(H, psi0, ti, tf, dt, print_nstep)
        use blas_wrappers
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        complex(dp), allocatable, intent(in) :: H(:,:)
        complex(dp), allocatable, intent(in) :: psi0(:)
        real(dp), intent(in) :: ti, tf, dt
        integer, intent(in) :: print_nstep
       
        complex(dp), allocatable :: psi_i(:)
        complex(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        complex(dp), allocatable, dimension(:) :: kt

        
        complex(dp) :: expt, norm, autocorr
        real(dp) :: t
        integer :: n, m
        integer :: tstep
        complex(dp) ::  alpha

        alpha=(0.0d0, -1.0d0)
        n = size(H, dim=2)
        m = size(psi0)
        ! some consistency checks
        if (n .ne. m) then
            print*, "!!ERROR: Dimensions of func and psi0 don't match."
            stop
        end if

        allocate(k1(n), k2(n), k3(n), k4(n))
        allocate(kt(n))
        allocate(psi_i(n))
    
        t = ti
        tstep = 0
        psi_i = psi0
        open(100, file='real_tprop.out')
        ! OUTFILE headers
        write(100,*) ""
        write(100,*) "   Real Time Propagation   "
        write(100,*) "---------------------------"
        write(100,*) ""
        
        ! writing for t0
        call zmul_zdotc(psi_i, psi_i, norm)
        ! expectation values
        call zmul_mv(H, psi_i, kt, alpha)
        call zmul_zdotc(psi_i,kt, expt)
        ! autocorr
        call zmul_zdotc(psi0, psi_i, autocorr)
        write(100,'(4f32.16)') t, abs(norm), abs(autocorr), abs(expt)  

        do while (t < tf)
            call zmul_mv(H, psi_i, k1, alpha)
            kt = (psi_i + (dt*0.5d0)*k1)
            call zmul_mv(H, kt, k2, alpha)
            kt = (psi_i + (dt*0.5d0)*k2)
            call zmul_mv(H, kt, k3, alpha)
            kt = (psi_i + (dt*1.0d0)*k3)
            call zmul_mv(H, kt, k4, alpha)
            psi_i = psi_i + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
            t = t + dt
            tstep = tstep + 1
            if (modulo(tstep, print_nstep) == 0) then
                ! norm
                call zmul_zdotc(psi_i, psi_i, norm)
                ! intermediate normalisation only for itp
                psi_i = 1/(norm)**(0.5d0) * psi_i
                call zmul_zdotc(psi_i, psi_i, norm)
                ! expectation values
                call zmul_mv(H, psi_i, kt, alpha)
                call zmul_zdotc(psi_i,kt, expt)
                ! autocorr
                call zmul_zdotc(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, abs(norm), abs(autocorr), abs(expt)  
            end if
        end do
        write(100,*) ""
        write(100,*) "   Summary of Real time propagation   "
        write(100,*) "--------------------------------------"
        write(100,'(a,i16)')   "ndim    = ", n
        write(100,'(a,f16.8)') "dt      = ", dt
        write(100,'(a,i16)')   "nstep   = ", tstep
        write(100,*) ""
        close(100)
    end subroutine real_tprop

    ! subroutine for imaginary time propagation
    subroutine imag_tprop(H, psi0, dt, print_nstep, etol)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        use blas_wrappers
        implicit none
        real(dp), allocatable, intent(in) :: H(:,:)
        real(dp), allocatable, intent(inout) :: psi0(:)
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
        real(dp) :: alpha
        
        n = size(H, dim=1)
        allocate(k1(n), k2(n), k3(n), k4(n))
        allocate(kt(n))
        allocate(psi_i(n))
        ! Input parameters
        t = 0.0d0
        tstep = 0
        psi_i = psi0
        E0 = 0.0d0
        call dmul_mv(H, psi_i, kt, alpha=1.0d0)
        call dmul_ddot(kt, psi_i, E0)
        Ei = E0
        open(100, file='imag_tprop.out')
        ! OUTFILE headers
        write(100,*) ""
        write(100,*) "Imaginary Time Propagation"
        write(100,*) "--------------------------"
        write(100,*) ""
        
        ! writing headers
        write(100,'(4a32)') 'Imag. time (a.u.)', '< psi_i | psi_i >', '< psi_i | psi_0 >', 'Energy (hartree)'
        do 
            ! print output to itp.out
            if (modulo(tstep, print_nstep) == 0) then
                ! norm
                call dmul_ddot(psi_i, psi_i, norm)
                ! autocorr
                call dmul_ddot(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, norm, autocorr, Ei 
            end if
            ! RK4 step
            call dmul_mv(H, psi_i, k1, alpha=-1.0d0)
            kt = (psi_i + (dt*0.5d0)*k1)
            call dmul_mv(H, kt, k2, alpha=-1.0d0)
            kt = (psi_i + (dt*0.5d0)*k2)
            call dmul_mv(H, kt, k3, alpha=-1.0d0)
            kt = (psi_i + (dt*1.0d0)*k3)
            call dmul_mv(H, kt, k4, alpha=-1.0d0)
            psi_i = psi_i + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
            t = t + dt
            tstep = tstep + 1
            ! intermediate normalisation only for itp
            call dmul_ddot(psi_i, psi_i, norm)
            psi_i = 1/(norm)**(0.5d0) * psi_i
            ! calculating energy
            call dmul_mv(H, psi_i, kt, alpha=1.0d0)
            call dmul_ddot(kt, psi_i, Ei)
            dE = abs(E0 - Ei)
            E0 = Ei
            if (dE < etol) then
                exit
            end if
        end do
        write(100,*) ""
        write(100,*) "Summary of Imaginary time propagation"
        write(100,*) "-------------------------------------"
        write(100,'(a,i16)')   "ndim    = ", n
        write(100,'(a,f16.8)') "dt      = ", dt
        write(100,'(a,i16)')   "nstep   = ", tstep
        write(100,'(a,f16.8)') "E_conv  = ", etol
        write(100,'(a,f16.8)') "Ei      = ", Ei
        write(100,*) ""
        close(100)

    end subroutine imag_tprop

end module rungekutta