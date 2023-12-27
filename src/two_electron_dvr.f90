module two_electron_dvr
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use omp_lib

    implicit none

    private

    public  :: print_te_dvr, oe_kinetic, oe_sc_single_well, oe_sc_multi_well, & 
               te_kinetic, te_sw_potential, te_mw_potential, &  
               oe_ho_hamiltonian, oe_scsw_hamiltonian, oe_scmw_hamiltonian, & 
               te_sw_hamiltonian, te_mw_hamiltonian, &
               sparse_te_sw_hamiltonian!, sparse_te_mw_hamiltonian

    ! grid and soft-coulomb potential parameters
    real(dp), public :: x0    ! grid size (-x0, x0)
    integer,  public :: n     ! number of grid points for each electron
    real(dp), public :: dx    ! grid spacing 
    real(dp), public :: m     ! mass of the electron
    real(dp), public :: z     ! atomic number?
    real(dp), public :: alpha ! single-well soft-coulomb parameter
    real(dp), public :: beta  ! e- correlation soft coloumb parameter
    real(dp), allocatable, public :: x(:)     ! grid points
    logical,  public :: multi_well_switch     ! set to true for multi-well
    logical,  public :: te_swtich             ! two-electron switch
    real(dp), allocatable, public :: atom_pos(:) ! position of nuclei
    real(dp), allocatable, public :: atom_num(:) ! atomic numbers
    
    real(dp), parameter, public :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 
    
    contains
    
    ! 1e- kinetic terms in sinc-dvr
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

    ! 1e- single-well soft-coulomb potential
    function oe_sc_single_well(i, j) &
        result(v_ij)
        implicit none
        integer,  intent(in) :: i,j      ! indices of 1e- matrix
        real(dp) :: v_ij                 ! matrix element V(i,j)
        if (i==j) then
            v_ij = -z/sqrt(x(i)**2 + alpha)
        else 
            v_ij = 0.0d0
        end if
    end function

    ! 1e- harmonic potential
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

    ! 1e- multi-well soft coulomb potential
    function oe_sc_multi_well(i, j) &
        result(v_ij)
        implicit none
        integer,  intent(in) :: i,j      ! indices of 1e- matrix
        real(dp) :: v_ij                 ! matrix element V(i,j)
        integer :: k, l
        l = size(atom_pos)
        if (i==j) then
            v_ij = 0.0d0
            do k=1,l
                v_ij = v_ij - atom_num(k)/sqrt((x(i)-atom_pos(k))**2 + alpha)
            end do
        else 
            v_ij = 0.0d0
        end if
    end function

    ! harmonic oscillator hamiltonian
    function oe_ho_hamiltonian(i,j) &
        result(h_ij)
        implicit none 
        integer, intent(in) :: i,j
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_harmonic_potential(i,j)
    end function oe_ho_hamiltonian

    ! 1e- hamiltonian for single-well soft-coulomb potential
    function oe_scsw_hamiltonian(i, j) &
        result(h_ij)
        integer, intent(in) :: i,j 
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_sc_single_well(i,j)
    end function oe_scsw_hamiltonian

    ! 1e- hamiltonian for multi-well soft-coulomb potential
    function oe_scmw_hamiltonian(i, j) &
        result(h_ij)
        integer, intent(in) :: i,j 
        real(dp) :: h_ij
        h_ij = oe_kinetic(i,j) + oe_sc_multi_well(i,j)
    end function oe_scmw_hamiltonian

    ! 2e- kinetic terms in sinc-dvr
    function te_kinetic(u, v) &
        result(t_uv)
        implicit none
        integer, intent(in) :: u,v ! 2e- matrix element indices
        real(dp) :: t_uv           ! matrix element T(u,v)
        integer :: i, j, k, l
        i = ((u-1)/n)+1
        j = ((v-1)/n)+1
        k = modulo((u-1),n)+1
        l = modulo((v-1),n)+1
        t_uv = 0.0d0
        if (k==l) then 
            t_uv = t_uv + oe_kinetic(i,j)
        end if 
        if (i==j) then
            t_uv = t_uv + oe_kinetic(k,l)
        end if
    end function te_kinetic

    ! 2e- single well soft-coulomb potential
    function te_sw_potential(u, v) &
        result(v_uv)
        implicit none
        integer,  intent(in) :: u,v             ! 2e- matrix element indices
        real(dp) :: v_uv                        ! matrix element V(u,v)
        integer :: i,j
        if (u == v) then
            i = ((u-1)/n)+1
            j = modulo((u-1),n)+1
            ! single electron terms
            v_uv = oe_sc_single_well(i, i) + oe_sc_single_well(j, j)
            ! electron-correlation term (note: has to be positive)
            v_uv = v_uv + (1.0d0)/sqrt(beta + (x(i)-x(j))**2)
        else
            v_uv = 0.0d0
        end if
    end function te_sw_potential

    ! 2e- multi well soft-coulomb potential
    function te_mw_potential(u, v) & 
        result(v_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: v_uv                        ! matrix element V(u,v)
        integer :: i,j
        if (u == v) then
            i = int((u-1)/n)+1
            j = modulo((u-1),n)+1
            ! single electron terms
            v_uv = oe_sc_multi_well(i, i) + oe_sc_multi_well(j, j)
            ! electron-correlation term (note: has to be positive)
            v_uv = v_uv + (1.0d0)/sqrt(beta + (x(i)-x(j))**2)
        else
            v_uv = 0.0d0
        end if
    end function te_mw_potential

    ! 2e- hamiltonian for single-well soft-coulomb potential
    function te_sw_hamiltonian(u, v) &
        result(h_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: h_uv                        ! matrix element H(u,v)
        h_uv = te_kinetic(u, v) + te_sw_potential(u, v)
    end function te_sw_hamiltonian

    ! 2e- hamiltonian for multi-well soft-coulomb potential
    function te_mw_hamiltonian(u, v) &
        result(h_uv)
        implicit none
        integer, intent(in) :: u,v              ! 2e- matrix element indices
        real(dp) :: h_uv                        ! matrix element H(u,v)
        h_uv = te_kinetic(u, v) + te_mw_potential(u, v)
    end function te_mw_hamiltonian

    ! sparse hamiltonian 
    subroutine sparse_te_sw_hamiltonian(h_array, row_idx, col_idx)
        implicit none
        real(dp), intent(inout) :: h_array(:)
        integer,  intent(out) :: row_idx(:)
        integer,  intent(out) :: col_idx(:)
        integer :: ndim, nsparse
        integer :: u, v
        integer :: i, j
        real(dp) :: atol
        real(dp) :: h_uv
        
        ndim = n**2
        atol = 1e-8
        nsparse= 2*n - 1
        ! allocate(h_array(ndim*nsparse))
        ! allocate(col_idx(ndim*nsparse))
        ! allocate(row_idx(ndim+1))
        !$omp parallel default(none) private(u,v,i,j,h_uv) shared(ndim, nsparse, atol, h_array, col_idx, row_idx)
        !$omp do
        do u=1,ndim
            i = 1
            row_idx(u) = (u-1)*nsparse + 1
            do v=1,ndim
                h_uv = te_sw_hamiltonian(u,v)
                if (abs(h_uv) >atol) then
                    j = (u-1)*nsparse + i
                    h_array(j) = h_uv
                    col_idx(j) = v
                    i = i + 1
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel
        row_idx(ndim+1) = ndim*nsparse + 1
        
    end subroutine sparse_te_sw_hamiltonian

    subroutine print_te_dvr(file_unit)
        integer, intent(in) :: file_unit 
        write(file_unit, *) ""
        write(file_unit, *) "    Two-Electron DVR    "
        write(file_unit, *) "------------------------"
        write(file_unit,'(a,f16.8)') " x0     = ", x0
        write(file_unit,'(a,f16.8)') " m      = ", m
        write(file_unit,'(a,f16.8)') " Z      = ", z
        write(file_unit,'(a,i16)')   " n      = ", n
        write(file_unit,'(a,f16.8)') " dx     = ", dx
        write(file_unit, *) "" 
        write(file_unit,'(a)') "Soft-Coulomb Parameters:"
        write(file_unit,'(a,f16.8)') " alpha  = ", alpha
        write(file_unit, '(a,f16.8)')" beta   = ", beta
        write(file_unit, *) "" 
    end subroutine print_te_dvr
    
end module two_electron_dvr

