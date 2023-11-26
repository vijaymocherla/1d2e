! ifort -c exact_2e_dvr.f90 -qopenmp 
module two_electron_dvr
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none

    private :: oe_kinetic, oe_sc_single_well, oe_sc_multi_well, &
               te_kinetic, te_sc_single_well, te_sc_multi_well

    public  :: oe_ho_hamiltonian, oe_scsw_hamiltonian, oe_scmw_hamiltonian, & 
               te_sw_hamiltonian, te_mw_hamiltonian, &
               imag_tprop, gen_trial_state, normalize, &
               omp_minus_ham_psi, omp_dotprod, omp_daxpy


    ! grid and soft-coulomb potential parameters
    real(dp), public :: x0    ! grid size (-x0, x0)
    integer,  public :: n     ! number of grid points for each electron
    integer,  public :: ndim  ! dimension of two-electron direct-product space 
    real(dp), public :: dx    ! grid spacing 
    real(dp), public :: m     ! mass of the electron
    real(dp), public :: alpha ! single-well soft-coulomb parameter
    real(dp), public :: beta  ! e- correlation soft coloumb parameter
    real(dp), allocatable, public :: x(:)     ! grid points
    logical,  public :: multi_well_switch     ! set to true for multi-well
    logical,  public :: te_swtich             ! two-electron switch
    real(dp), allocatable, public :: wells(:) ! position of nuclei
    
    real(dp), parameter, public :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 
    
    contains
    
    function oe_kinetic(i, j) &
        result(t_ij)
        implicit none
        integer, intent(in) :: i,j ! indices of 1e- matrix
        real(dp) :: t_ij           ! matrix element T(i,j)
        ! matrix elements
        if (i==j) then
            t_ij =  (-1)**(i-i) * (1.0d0 / (2.0d0 * m * dx**2)) * pi_dp**2 / 3.0d0
        else
            t_ij = (-1)**(i-j) * (1.0d0 / (2.0d0 * m * dx**2)) * 2.0d0 / (i-j)**2
        end if
    end function oe_kinetic

    function te_kinetic(u, v) &
        result(t_uv)
        implicit none
        integer, intent(in) :: u,v ! 2e- matrix element indices
        real(dp) :: t_uv           ! matrix element T(u,v)
        t_uv = oe_kinetic(int(u/n)+1,int(v/n)+1) * oe_kinetic(modulo(u,n)+1,modulo(v,n)+1) 
    end function te_kinetic

    function oe_sc_single_well(i, j) &
        result(v_ij)
        implicit none
        integer,  intent(in) :: i,j      ! indices of 1e- matrix
        real(dp) :: v_ij                 ! matrix element V(i,j)
        if (i==j) then
            v_ij = -1.0d0/sqrt(x(i)**2 + alpha)
        else 
            v_ij = 0.0d0
        end if
    end function

    function oe_harmonic_potential(i, j) &
        result(v_ij)
        implicit none
        integer,  intent(in) :: i,j      ! indices of 1e- matrix
        real(dp) :: v_ij                 ! matrix element V(i,j)
        if (i==j) then
            v_ij = 0.5d0 * x(i)**2
        else 
            v_ij = 0.0d0
        end if
    end function

    function oe_sc_multi_well(i, j) &
        result(v_ij)
        implicit none
        integer,  intent(in) :: i,j      ! indices of 1e- matrix
        real(dp) :: v_ij                 ! matrix element V(i,j)
        integer :: k, n
        n = size(wells)
        if (i==j) then
            v_ij = 0.0d0
            do k=1,n
                v_ij = v_ij - 1.0d0/sqrt((x(i)-wells(k))**2 + alpha)
            end do
        else 
            v_ij = 0.0d0
        end if
    end function

    function oe_ho_hamiltonian(i,j) &
        result(h_ij)
        implicit none 
        integer, intent(in) :: i,j
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_harmonic_potential(i,j)
    end function oe_ho_hamiltonian

    function oe_scsw_hamiltonian(i, j) &
        result(h_ij)
        integer, intent(in) :: i,j 
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_sc_single_well(i,j)
    end function oe_scsw_hamiltonian

    function oe_scmw_hamiltonian(i, j) &
        result(h_ij)
        integer, intent(in) :: i,j 
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_sc_multi_well(i,j)
    end function oe_scmw_hamiltonian

    function te_sc_single_well(u, v) &
        result(v_uv)
        implicit none
        integer,  intent(in) :: u,v             ! 2e- matrix element indices
        real(dp) :: v_uv                        ! matrix element V(u,v)
        integer :: i,j
        if (u == v) then
            i = int(u/n) + 1
            j = modulo(u, n) + 1
            ! single electron terms
            v_uv = oe_sc_single_well(i, i) * oe_sc_single_well(j, j)
            ! electron-correlation term
            v_uv = v_uv + (-1.0d0)/sqrt(beta + (x(i)-x(j))**2)
        else
            v_uv = 0.0d0
        end if
    end function te_sc_single_well

    function te_sc_multi_well(u, v) & 
        result(v_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: v_uv                        ! matrix element V(u,v)
        integer :: i,j
        if (u == v) then
            i = int(u/n) + 1
            j = modulo(u, n) + 1
            ! single electron terms
            v_uv = oe_sc_multi_well(i, i) * oe_sc_multi_well(j, j)
            ! electron-correlation term
            v_uv = v_uv + (-1.0d0)/sqrt(beta + (x(i)-x(j))**2)
        else
            v_uv = 0.0d0
        end if
    end function te_sc_multi_well
    
    function te_sw_hamiltonian(u, v) &
        result(h_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: h_uv                        ! matrix element H(u,v)
        h_uv = te_kinetic(u, v) + te_sc_single_well(u, v)
    end function te_sw_hamiltonian

    function te_mw_hamiltonian(u, v) &
        result(h_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: h_uv                        ! matrix element H(u,v)
        h_uv = te_kinetic(u, v) + te_sc_multi_well(u, v)
    end function te_mw_hamiltonian


    subroutine omp_minus_ham_psi(a, psi_i, psi_j)
        use omp_lib
        implicit none
        real(dp), intent(in) :: a
        real(dp), dimension(:), intent(in)  :: psi_i
        real(dp), dimension(:), intent(out) :: psi_j
        integer :: u, v
        real(dp) :: tmp
        psi_j = 0.0d0
        if (te_swtich) then
            if (multi_well_switch) then
                !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
                !$omp do
                do u=1,ndim
                    tmp = 0.0d0
                    do v=1,ndim
                        tmp = tmp + a*te_mw_hamiltonian(u,v)*psi_i(v)
                    end do
                    psi_j(u) = tmp
                end do
                !$omp end do
                !$omp end parallel
            else
                !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
                !$omp do 
                do u=1,ndim
                    tmp = 0.0d0
                    do v=1,ndim
                        tmp = tmp + a*te_sw_hamiltonian(u,v)*psi_i(v)
                    end do
                    psi_j(u) = tmp
                end do
                !$omp end do
                !$omp end parallel
            end if
        else
            !$omp parallel shared(ndim, psi_i, psi_j) private(u,v)
            !$omp do
            do u=1,ndim
                tmp = 0.0d0
                do v=1,ndim
                    tmp = tmp + a*oe_ho_hamiltonian(u,v)*psi_i(v)
                end do
                psi_j(u) = tmp
            end do
            !$omp end do
            !$omp end parallel
        end if
    end subroutine


    subroutine omp_dotprod(psi_i, psi_j, result)
        use omp_lib
        implicit none
        real(dp), dimension(:), intent(in)  :: psi_i, psi_j
        real(dp), intent(out) :: result
        integer :: u
        result = 0.0d0
        !$omp parallel shared(ndim, psi_i, psi_j) private(u)
        !$omp do reduction ( + : result )
        do u = 1, ndim
            result = result + psi_i(u) * psi_j(u)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine


    subroutine omp_daxpy(a, x, y, z)
        use omp_lib
        implicit none
        real(dp), intent(in) :: a
        real(dp), dimension(:), intent(in)  :: x, y
        real(dp), dimension(:), intent(out) :: z
        integer :: i
        z = 0.0d0
        !$omp parallel shared(ndim, a, x, y) private(i)
        !$omp do
        do i=1,ndim
            z(i) = a*x(i) + y(i)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine


    subroutine normalize(psi)
        real(dp), intent(inout) :: psi(:)
        integer :: u
        real(dp) :: norm, inv_norm
        call omp_dotprod(psi, psi, norm)
        inv_norm = 1.0d0/(norm)**(0.5)
        !$omp parallel shared(ndim, inv_norm, psi) private(u)
        !$omp do
        do u=1,ndim
            psi(u) = inv_norm * psi(u)
        end do
        !$omp end do 
        !$omp end parallel
    end subroutine normalize


    subroutine gen_trial_state(psi)
        real(dp), intent(out) :: psi(:)
        integer :: u, i, j
        real(dp) :: zeta, rij
        zeta = 0.5 ! effective nuclear charge for helium in 3D
        ! generating trail state psi(x1,x2) = |x1-x2| exp(-x1)*exp(-x2)
        !$omp parallel shared(ndim, x, psi) private(u)
        !$omp do
        do u=1,ndim
            i = int(u/n) + 1
            j = modulo(u, n) + 1
            rij = x(i)-x(j)
            psi(u) = rij*exp(-zeta*abs(x(i)))*exp(-zeta*abs(x(j)))
        end do
        !$omp end do
        !$omp end parallel
        ! normalizing trial state
        call normalize(psi) 
    end subroutine


    subroutine imag_tprop(psi0, dt, print_nstep, etol, Ei, tstep)
        implicit none
        real(dp), allocatable, intent(in) :: psi0(:)
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
        integer :: ndim
        real(dp) :: a

        ! allocating arrays
        ndim = n**2
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
            call normalize(psi_i)
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

end module two_electron_dvr

