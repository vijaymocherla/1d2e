module itprop
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    use blas_wrappers, only: csr_dmul_mv
    use helpers, only: omp_normalize, omp_daxpy, omp_dotprod
    use two_electron_dvr

    implicit none
    private 

    public :: itp_on_the_fly, itp_sparse, gen_trial_state

    contains
    
    ! generates a trial wave fxn for ITP
    subroutine gen_trial_state(psi)
        real(dp), intent(out) :: psi(:)
        integer :: u, i, j
        real(dp) :: zeta
        integer :: ndim
        zeta = 1.0 ! effective nuclear charge for helium in 3D
        ndim = size(psi)
        ! generating trail state psi(x1,x2) = |x1-x2| exp(-x1)*exp(-x2)
        !$omp parallel shared(ndim, x, psi) private(u)
        !$omp do
        do u=1,ndim
            i = ((u-1)/n) + 1
            j = modulo(u-1, n) + 1
            psi(u) = exp(-zeta*(abs(x(i))+ abs(x(j))))
        end do
        !$omp end do
        !$omp end parallel
        ! normalizing trial state
        call omp_normalize(psi) 
    end subroutine
    
    ! parallelised routine for psi_j = a * hamiltonian * psi_i
    subroutine omp_minus_ham_psi(a, psi_i, psi_j)
        implicit none
        real(dp), intent(in) :: a
        real(dp), dimension(:), intent(in)  :: psi_i
        real(dp), dimension(:), intent(out) :: psi_j
        integer :: u, v, ndim
        
        ndim = size(psi_i)
        psi_j = 0.0d0
        if (te_swtich) then
            if (multi_well_switch) then
                !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
                !$omp do
                do u=1,ndim
                    do v=1,ndim
                        psi_j(u) = psi_j(u) + a*te_mw_hamiltonian(u,v)*psi_i(v)
                    end do
                end do
                !$omp end do
                !$omp end parallel
            else
                !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
                !$omp do 
                do u=1,ndim
                    do v=1,ndim
                        psi_j(u) = psi_j(u) + a*te_sw_hamiltonian(u,v)*psi_i(v)
                    end do
                end do
                !$omp end do
                !$omp end parallel
            end if
        else
            !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
            !$omp do
            do u=1,ndim
                do v=1,ndim
                    psi_j(u) = psi_j(u) + a*oe_scsw_hamiltonian(u,v)*psi_i(v)
                end do
            end do
            !$omp end do
            !$omp end parallel
        end if
    end subroutine

    ! imaginary time propagation for finding the ground state wave fxn  
    ! with on-the-fly calculations of hamiltonian elements 
    subroutine itp_on_the_fly(psi0, dt, print_nstep, etol, Ei, tstep)
        implicit none
        real(dp), allocatable, intent(inout) :: psi0(:)
        integer,  intent(in)    :: print_nstep
        integer, intent(inout)  :: tstep
        real(dp), intent(inout) :: Ei
        real(dp), intent(in)    :: dt
        real(dp), intent(in)    :: etol
        real(dp) :: t  ! time
        real(dp), allocatable, dimension(:) :: psi_i
        real(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        real(dp), allocatable, dimension(:) :: kt
        real(dp) :: const ! constant for daxpy operations
        real(dp) :: norm, autocorr, E0, dE
        real(dp) :: a
        integer :: ndim
        ! allocating arrays
        ndim = size(psi0)
        a = -1.0d0
        allocate(k1(ndim), k2(ndim), k3(ndim), k4(ndim))
        allocate(kt(ndim))
        allocate(psi_i(ndim))
        
        ! initiating some variables
        t = 0.0d0
        tstep = 0
        const = 0.0d0
        psi_i = psi0
        E0 = 0.0d0
        call omp_minus_ham_psi(a, psi_i, kt) ! omp_minus_ham_psi gives -(\hat{H}.\psi)
        call omp_dotprod(kt, psi_i, Ei)   ! Ei = <\psi|-\hat{H}|\psi>
        Ei = -1.0d0*Ei                    ! multiply -1.0d0 to get correct energy
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
                call omp_dotprod(psi_i, psi_i, norm)
                ! autocorr
                call omp_dotprod(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, norm, autocorr, Ei 
            end if
            
            ! RK4 step
            call omp_minus_ham_psi(a, psi_i, k1)
            call omp_daxpy(0.5d0*dt, k1, psi_i, kt)
            call omp_minus_ham_psi(a, kt, k2)
            call omp_daxpy(0.5d0*dt, k2, psi_i, kt)
            call omp_minus_ham_psi(a, kt, k3)
            call omp_daxpy(1.0d0*dt, k3, psi_i, kt)
            call omp_minus_ham_psi(a, kt, k4)
            ! parallelising: psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            call omp_daxpy(dt/6.0d0, k1, psi_i, kt)
            call omp_daxpy(dt/3.0d0, k2, kt, psi_i)
            call omp_daxpy(dt/3.0d0, k3, psi_i, kt)
            call omp_daxpy(dt/6.0d0, k4, kt, psi_i)
            ! intermediate normalisation only for itp
            call omp_normalize(psi_i)
            ! time step increment
            t = t + dt
            tstep = tstep + 1
            ! calculating energy (Ei)
            call omp_minus_ham_psi(a, psi_i, kt) ! omp_minus_ham_psi gives -(\hat{H}.\psi)
            call omp_dotprod(kt, psi_i, Ei)   ! Ei = <\psi|-\hat{H}|\psi>
            Ei = -1.0d0*Ei                    ! multiply -1.0d0 to get correct energy
            dE = abs(E0 - Ei)
            E0 = Ei
            ! checking for Energy convergence
            if (dE < etol) then
                ! norm
                call omp_dotprod(psi_i, psi_i, norm)
                ! autocorr
                call omp_dotprod(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, norm, autocorr, Ei
                psi0 = psi_i
                exit
            end if
        end do    
        write(100,*) ""
        write(100,*) "Summary of Imaginary time propagation"
        write(100,*) "-------------------------------------"
        write(100,'(a,i16)')   "n       = ", n
        write(100,'(a,f16.8)') "dt      = ", dt
        write(100,'(a,i16)')   "nstep   = ", tstep
        write(100,'(a,f16.8)') "E_conv  = ", etol
        write(100,'(a,f16.8)') "Ei      = ", Ei
        write(100,*) ""
        close(100)
    end subroutine itp_on_the_fly
    

    ! imaginary time propagation for a sparse hamiltonian
    subroutine itp_sparse(h_array, h_row, h_col, psi0, dt, print_nstep, etol, Ei, tstep)
        implicit none
        real(dp), allocatable, intent(inout) :: h_array(:)
        integer(dp),  allocatable, intent(in) :: h_row(:)
        integer(dp),  allocatable, intent(in) :: h_col(:)
        real(dp), allocatable, intent(inout) :: psi0(:)
        integer,  intent(in)    :: print_nstep
        integer, intent(inout)  :: tstep
        real(dp), intent(inout) :: Ei
        real(dp), intent(in)    :: dt
        real(dp), intent(in)    :: etol
        real(dp) :: t  ! time
        real(dp), allocatable, dimension(:) :: psi_i
        real(dp), allocatable, dimension(:) :: k1, k2, k3, k4 
        real(dp), allocatable, dimension(:) :: kt
        real(dp) :: const ! constant for daxpy operations
        real(dp) :: norm, autocorr, E0, dE
        real(dp) :: a
        integer :: ndim

        ! allocating arrays
        ndim = size(psi0)
        a = -1.0d0
        allocate(k1(ndim), k2(ndim), k3(ndim), k4(ndim))
        allocate(kt(ndim))
        allocate(psi_i(ndim))
        
        ! initiating some variables
        t = 0.0d0
        tstep = 0
        const = 0.0d0
        psi_i = psi0
        E0 = 0.0d0
        ! changing H = -1.0 * H
        h_array = -1.0d0* h_array
        ! omp_minus_ham_psi gives -(\hat{H}.\psi)
        call csr_dmul_mv(h_array, h_row, h_col, psi_i, kt)
        call omp_dotprod(kt, psi_i, Ei)   ! Ei = <\psi|-\hat{H}|\psi>
        Ei = -1.0d0*Ei                    ! multiply -1.0d0 to get correct energy
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
                call omp_dotprod(psi_i, psi_i, norm)
                ! autocorr
                call omp_dotprod(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, norm, autocorr, Ei 
            end if
            
            ! RK4 step
            ! computing k1
            call csr_dmul_mv(h_array, h_row, h_col, psi_i, k1)
            ! computing k2
            call omp_daxpy(0.5d0*dt, k1, psi_i, kt)
            call csr_dmul_mv(h_array, h_row, h_col, kt, k2)
            ! computing k3
            call omp_daxpy(0.5d0*dt, k2, psi_i, kt)
            call csr_dmul_mv(h_array, h_row, h_col, kt, k3)
            ! computing k4
            call omp_daxpy(1.0d0*dt, k3, psi_i, kt)
            call csr_dmul_mv(h_array, h_row, h_col, kt, k4)
            ! parallelising: psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            call omp_daxpy(dt/6.0d0, k1, psi_i, kt)
            call omp_daxpy(dt/3.0d0, k2, kt, psi_i)
            call omp_daxpy(dt/3.0d0, k3, psi_i, kt)
            call omp_daxpy(dt/6.0d0, k4, kt, psi_i)
            ! intermediate normalisation only for itp
            call omp_normalize(psi_i)
            ! time step increment
            t = t + dt
            tstep = tstep + 1
            
            ! calculating energy (Ei)
            call csr_dmul_mv(h_array, h_row, h_col, psi_i, kt)
            ! call omp_minus_ham_psi(a, psi_i, kt) ! omp_minus_ham_psi gives -(\hat{H}.\psi)
            call omp_dotprod(kt, psi_i, Ei)   ! Ei = <\psi|-\hat{H}|\psi>
            Ei = -1.0d0*Ei                    ! multiply -1.0d0 to get correct energy
            dE = abs(E0 - Ei)
            E0 = Ei
            ! checking for Energy convergence
            if (dE < etol) then
                ! norm
                call omp_dotprod(psi_i, psi_i, norm)
                ! autocorr
                call omp_dotprod(psi0, psi_i, autocorr)
                write(100,'(4f32.16)') t, norm, autocorr, Ei
                psi0 = psi_i
                exit
            end if
        end do    
        write(100,*) ""
        write(100,*) "Summary of Imaginary time propagation"
        write(100,*) "-------------------------------------"
        write(100,'(a,i16)')   "n       = ", n
        write(100,'(a,f16.8)') "dt      = ", dt
        write(100,'(a,i16)')   "nstep   = ", tstep
        write(100,'(a,f16.8)') "E_conv  = ", etol
        write(100,'(a,f16.8)') "Ei      = ", Ei
        write(100,*) ""
        close(100)
    end subroutine itp_sparse
end module itprop