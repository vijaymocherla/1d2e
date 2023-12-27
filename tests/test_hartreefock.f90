program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use lapack_wrappers, only:eigsh
    use hartreefock

    implicit none
    real(dp), allocatable :: h_core(:,:)
    real(dp), allocatable :: v_ee(:,:)  
    real(dp), allocatable :: epsilon(:)
    real(dp), allocatable :: mo_coeff(:,:)
    integer :: i, j, print_nstep
    integer :: maxiter
    real(dp) :: etol, Ei
    integer :: ci
    character (len=32) :: arg
    
    ci = 1
    ! Reading input parameters
    n = 8                        ! size of 1e- grids
    x0 = 10.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    etol = 0.00000000001         ! energy absolute tolerance 
    maxiter = 200                ! max scf cycle
    n_el = 2                     ! no. of electrons
    damping_factor = 0.0d0       ! damping factor  
    z = 2.0d0

    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-n") then
            call get_command_argument(ci+1,arg)
            read(arg, '(i32)') n
            ci = ci + 2
        else if (trim(arg)=="-x0") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') x0
            ci = ci + 2
        else if (trim(arg)=="-z") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') z
            ci = ci + 2
        else if (trim(arg)=="-damping_factor") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') damping_factor
            ci = ci + 2
        else if (trim(arg)=="-alpha") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') alpha
            ci = ci + 2
        else if (trim(arg)=="-beta") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') beta
            ci = ci + 2
        else if (trim(arg)=="-etol") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') etol
            ci = ci + 2    
        else if (trim(arg)=="-print_nstep") then
            call get_command_argument(ci+1,arg)
            read(arg, '(i32)') print_nstep
            ci = ci + 2
        else
            exit
        end if        
    end do
    print*, z
    ! setting parameters for grids and calculations
    m = 1.0d0                    ! mass of the electron
    dx = 2.0d0*x0 / real(n-1,8)  ! grid-spacing
    multi_well_switch = .false.  ! default for single well
    te_swtich = .true.           ! switch to test one-electron hamiltonian 
    
    print*, ""
    print*, "    Two-Electron DVR    "
    print*, "------------------------"
    print('(a,f16.8)'), "x0     = ", x0
    print('(a,f16.8)'), "m      = ", m
    print('(a,f16.8)'), "Z      = ", z
    print('(a,i16)'),   "n      = ", n
    print('(a,f16.8)'), "dx     = ", dx
    print*, "" 
    print('(a)'), "Soft-Coulomb Parameters:"
    print('(a,f16.8)'), "alpha  = ", alpha
    print('(a,f16.8)'), "beta   = ", beta
    print*, "" 
    ! Position grid points
    allocate(x(n))
    do i=1,n
        x(i) = -x0 + (i-1)*dx
    end do
    ! Allocating and initiating array
    allocate(h_core(n,n), v_ee(n,n))
    allocate(mo_coeff(n,n), epsilon(n))
    h_core = 0.0d0
    v_ee = 0.0d0
    print*, "Generating Hamiltonian Matrix......"
    print*, "" 
    ! core hamiltonian matrix
    do i=1,n
        h_core(i,i) = oe_kinetic(i,i) + oe_sc_single_well(i,i)
        do j=i+1,n
            h_core(i,j) = oe_kinetic(i,j)
            h_core(j,i) = h_core(i,j)
        end do
    end do
    ! electron-electron correlation
    do i=1,n
        v_ee(i,i) = (1.0d0)/sqrt(beta)
        do j=i+1,n
            v_ee(i,j) = 1.0d0/sqrt(beta + (x(i)-x(j))**2)
            v_ee(j,i) = v_ee(i,j)
        end do
    end do
    print*, "Running SCF........"
    ! run scf
    call scf_cycle(h_core, v_ee, etol, maxiter, Ei, mo_coeff, epsilon)
    print*, "Completed running SCF!!!"
    


end program main