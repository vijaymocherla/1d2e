module real_tprop
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use blas_wrappers, only: csr_zmul_mv, zmul_zdotc, dmul_ddot
    use helpers, only: omp_normalize, omp_zaxpy, omp_dotprod
    use two_electron_dvr
    implicit none

    private

    public :: rtp_sparse, gen_intial_wfn
    
    contains 

    subroutine gen_intial_wfn(psi, p01, p02, zeta)
        real(dp),  intent(in) :: p01, p02, zeta
        complex(dp), intent(out) :: psi(:)
        integer :: u, i, j
        real(dp) :: norm
        complex(dp) :: res
        complex(dp) :: c1, c2
        integer :: ndim
        ndim = size(psi)
        c1 = cmplx(zeta, -p01, kind=dp)
        c2 = cmplx(zeta, -p02, kind=dp)
        ! generating trail state
        !$omp parallel shared(ndim, x, psi, c1, c2) private(u)
        !$omp do
        do u=1,ndim
            i = ((u-1)/n) + 1
            j = modulo(u-1, n) + 1
            psi(u) = exp(-c1*abs(x(i)) - c2*abs(x(j)))
        end do
        !$omp end do
        !$omp end parallel
        ! normalizing trial state
        call zmul_zdotc(psi, psi, res)
        norm = real(res,kind=dp)
        psi = 1.0d0/(norm)**(0.5) * psi
    end subroutine gen_intial_wfn

    subroutine rtp_sparse(h_array, h_row, h_col, psi0, ti, tf, dt, print_nstep, Ei, tstep)
        implicit none
        complex(dp), allocatable, intent(inout) :: h_array(:)
        integer, allocatable, intent(in) :: h_row(:)
        integer, allocatable, intent(in) :: h_col(:)
        complex(dp), allocatable, intent(inout) :: psi0(:)
        real(dp),  intent(in)      :: ti, tf
        integer,   intent(in)      :: print_nstep
        integer,   intent(inout)   :: tstep
        real(dp), intent(inout)    :: Ei
        complex(dp) :: res
        real(dp), intent(in) :: dt
        real(dp) :: norm, autocorr
        real(dp) :: t  ! time
        complex(dp), allocatable, dimension(:) :: kt
        complex(dp), allocatable, dimension(:) :: psi_i
        complex(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        complex(dp) :: const ! constant for zaxpy operations
        character (len=32) :: wfn_filename
        integer :: ndim, i, k
        ! allocating arrays
        ndim = size(psi0)
        allocate(k1(ndim), k2(ndim), k3(ndim), k4(ndim))
        allocate(kt(ndim))

        allocate(psi_i(ndim))
        
        ! initiating some variables
        h_array  = cmplx(0.0d0,-1.0d0, kind=dp)*h_array
        t = ti
        tstep = 0
        k = 0 
        psi_i = psi0
        open(100, file='real_tprop.out')
        ! OUTFILE headers
        write(100,*) ""
        write(100,*) "   Real Time Propagation   "
        write(100,*) "---------------------------"
        write(100,*) ""
        ! writing headers
        write(100,'(4a32)') 'time (a.u.)', '< psi_i | psi_i >', '< psi_i | psi_0 >', 'Energy (hartree)'
        do while(t < tf)
            ! print output to itp.out
            if (modulo(tstep, print_nstep) == 0) then
                ! norm
                call zmul_zdotc(psi_i, psi_i, res)
                norm = real(res, kind=dp)
                ! autocorr
                call zmul_zdotc(psi0, psi_i, res)
                autocorr = abs(res)**2
                ! calculating energy Ei = <\psi|-\hat{H}|\psi>
                call csr_zmul_mv(h_array, h_row, h_col, psi_i, kt)
                call zmul_zdotc(kt, psi_i, res)
                Ei = imag(res)

                write(100,'(4f32.16)') t, norm, autocorr, Ei 
                write(wfn_filename,'(1a3,1i3.3,1a4)') 'rtp',k,'.wfn'
                open(200, file=trim(wfn_filename))
                    do i=1,ndim
                        write(200,*) psi_i(i) 
                    end do    
                close(200)
                k = k + 1
            end if
            ! RK4 step
            ! computing k1
            call csr_zmul_mv(h_array, h_row, h_col, psi_i, k1)
            ! computing k2
            const = cmplx(0.5d0*dt, 0.0d0, kind=dp)
            call omp_zaxpy(const, k1, psi_i, kt)
            call csr_zmul_mv(h_array, h_row, h_col, kt, k2)
            ! computing k3
            call omp_zaxpy(const, k2, psi_i, kt)
            call csr_zmul_mv(h_array, h_row, h_col, kt, k3)
            ! computing k4
            const = cmplx(1.0d0*dt, 0.0d0, kind=dp)
            call omp_zaxpy(const, k3, psi_i, kt)
            call csr_zmul_mv(h_array, h_row, h_col, kt, k4)
            ! parallelising: psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            !$omp parallel private(i) shared(ndim, dt, psi_i, k1, k2, k3, k4)
            !$omp do
            do i=1, ndim
                psi_i(i) = psi_i(i) + dt/6.0d0*k1(i) + dt/3.0d0*k2(i) + dt/3.0d0*k3(i) + dt/6.0d0*k4(i)
            end do
            !$omp end do
            !$omp end parallel

            ! time step increment
            t = t + dt
            tstep = tstep + 1
            
        end do 
        psi0 = psi_i   
        ! write(100,*) ""
        ! write(100,*) "   Summary of Real time propagation   "
        ! write(100,*) "--------------------------------------"
        ! write(100,'(a,i16)')   "n               = ", n
        ! write(100,'(a,i16)')   "print_nstep     = ", print_nstep
        ! write(100,'(a,f16.8)') "dt              = ", dt
        ! write(100,'(a,f16.8)') "ti              = ", ti
        ! write(100,'(a,f16.8)') "tf              = ", tf
        ! write(100,'(a,f16.8)') "Ei              = ", Ei
        ! write(100,*) ""
        close(100)
    end subroutine rtp_sparse


end module real_tprop