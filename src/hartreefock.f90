! hartreefock.f90
!
! TO DO:
! - Change guess density    
! - parallelise elementwise multiplication in comp_two_electron_integral
! - Move to checking density for convergence
! - Check if DIIS is needed?
!
!
module hartreefock
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use omp_lib
    use lapack_wrappers, only:eigsh
    use blas_wrappers, only: dmul_mm
    use two_electron_dvr, only: print_te_dvr
    implicit none

    private :: comp_two_electron_integral, calc_energy
    
    public :: scf_cycle

    logical, public  :: scf_damping_switch
    logical, public  :: scf_diis_switch
    real(dp), public :: damping_factor
    integer, public  :: n_el 

    contains
    
    subroutine comp_two_electron_integral(v_ee, p, g)
        real(dp), allocatable, intent(in) :: v_ee(:,:)
        real(dp), allocatable, intent(in) :: p(:,:)
        real(dp), allocatable, intent(out):: g(:,:)
        integer :: n
        integer :: u, v 
        real(dp):: const
        const = -0.50d0
        n = size(v_ee, dim=1)
        g = -0.5d0*p*v_ee ! exchange term

        !$omp parallel private(u) shared(n, g, p, v_ee)
        !$omp do
        do u=1,n
            do v=1,n
                g(u,u) = g(u,u) + v_ee(u,v)*p(v,v)  ! coulomb term
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine comp_two_electron_integral

    subroutine calc_energy(h_core, f, p, Ei)
        implicit none
        real(dp), allocatable, intent(in) :: h_core(:,:)
        real(dp), allocatable, intent(in) :: f(:,:)
        real(dp), allocatable, intent(in) :: p(:,:)
        real(dp), intent(out) :: Ei
        real(dp) :: en, tmp
        integer :: n
        integer :: u,v
        n = size(h_core, dim=1)
        en = 0.0d0
        !$omp parallel private(u) shared(n, h_core, f)
        !$omp do reduction( + : en)
        do u=1,n
            tmp = 0.0d0
            do v=1,n
                tmp = tmp + 0.5d0*(h_core(u,v) + f(u,v))*p(v,u)
            end do
            en = en + tmp
        end do
       !$omp end do
       !$omp end parallel
        Ei = en
    end subroutine calc_energy
    
    subroutine scf_cycle(h_core, v_ee, etol, maxiter, Ei, mo_coeff, epsilon)
        implicit none
        real(dp), allocatable, intent(in) :: h_core(:,:)
        real(dp), allocatable, intent(in) :: v_ee(:,:)
        real(dp), allocatable, intent(inout) :: mo_coeff(:,:)
        real(dp), allocatable, intent(inout) :: epsilon(:)
        real(dp), intent(inout) :: Ei
        real(dp), intent(in) :: etol
        integer,  intent(in) :: maxiter

        real(dp), allocatable :: g(:,:)
        real(dp), allocatable :: p(:,:)
        real(dp), allocatable :: p_tmp(:,:), p_new(:,:)
        real(dp), allocatable :: f(:,:)
        real(dp), allocatable :: occ(:,:)
        integer :: i, n
        real(dp) :: dE, E0
        logical :: convergence

        n = size(h_core, dim=1)
        open(100, file='scf_cycle.out')
        call print_te_dvr(100)
        write(100, *) ""
        write(100, *) "       Hartree-Fock       "
        write(100, *) "--------------------------"
        write(100,'(a,i16)') " no. of orbitals        = ", n
        write(100,'(a,i16)') " no. of electrons       = ", n_el
        write(100,'(a,l8)')  " scf_damping_switch     = ", scf_damping_switch
        write(100,'(a,l8)')  " scf_diis_switch        = ", scf_diis_switch
        write(100,'(a,f12.8)') " damping_factor         = ", damping_factor
        write(100,'(a,f18.16)')" e_tolerance            = ", etol
        write(100,'(a,i16)') " maximum scf-iterations = ", maxiter

        write(100, *) ""         
        
        ! allocating and initiating arrays
        allocate(f(n,n), g(n,n), p(n,n), p_tmp(n,n), p_new(n,n), occ(n,n))
        ! orbital occupation number
        occ = 0.0d0
        do i=1,n_el/2
            occ(i,i) = 2.0d0
        end do
        write(100,*) "                      SCF iterations                      "
        write(100,*) "----------------------------------------------------------" 
        write(100,*) ""
        write(100,'(a8,a20,a20)') "S.No", "Energy (hartree)", "dE (hartree)"
        write(100,*) "----------------------------------------------------------" 
        Ei = 0.0d0
        ! guess density
        p = 0.0d0
        g = 0.0d0
        p_new = 0.0d0
        p_tmp = 0.0d0
        convergence = .false.
        E0 = 1.0d0
        ! SCF iteration cycle
        do i=1,maxiter
            ! compute the two-electron term for the new density
            call comp_two_electron_integral(v_ee, p, g)
            f = h_core + g
            ! calculate energy
            call calc_energy(h_core, f, p, Ei)
            write(100,'(i8,f20.12,f20.12)') i, Ei, dE 
            ! check for convergence
            dE = abs(Ei - E0)
            if (dE < etol) then
                convergence = .true.
                exit
            end if
            E0  = Ei
            ! diagonalizing the fock-matrix
            call eigsh(f, epsilon, mo_coeff)
            ! calculating the new density
            call dmul_mm(occ, transpose(mo_coeff),p_tmp)
            call dmul_mm(mo_coeff, p_tmp, p_new)
            ! damping the density
            if (scf_damping_switch) then
                p = (1.0d0 - damping_factor)*p_new + damping_factor * p
            else
                p = p_new
            end if
        end do
        deallocate(g, p, p_tmp, p_new, occ)
        write(100,*) "----------------------------------------------------------" 
        write(100,*) " "
        if (convergence) then
            write(100,*) ""
            write(100,*) ""
            write(100,*) "---------------------------------------"
            write(100,*) "    Molecular Orbital Energies (Eh)    "
            write(100,*) "---------------------------------------"
            do i=1,n
                write(100,'(i8,f16.8)') i, epsilon(i)
            end do
            write(100,*) "---------------------------------------"
            write(100,*) ""
            write(100,'(a,f32.16)') "SCF converged! Final Ground-State Energy (Eh): ", Ei
        else 
            write(100,*) "!!! ERROR: SCF did NOT converge, maximum no. of iterations exceeded." 
        end if
        close(100)
    end subroutine scf_cycle
    
    

end module hartreefock