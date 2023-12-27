!
!   Imaginary time propagation of the wavefunction 
!
program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use lapack_wrappers, only:eigsh
    use imag_tprop
    use helpers, only:write_wfn

    implicit none
    real(dp), allocatable :: h_array(:) 
    integer,  allocatable :: h_row(:)
    integer,  allocatable :: h_col(:)
    real(dp), allocatable :: psi0(:)
    complex(dp), allocatable :: psi(:)
    integer  :: i, k, tstep, print_nstep, ndim, nsparse
    real(dp) :: dt
    real(dp) :: etol, Ei
    integer  :: ci
    character (len=32) :: arg
    character (len=128) :: comment

    ci = 1
    ! Reading input parameters
    n = 32                       ! size of 1e- grids
    x0 = 15.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    etol = 0.000000001           ! energy absolute tolerance 
    dt  = 0.001d0
    print_nstep = 100
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
        else if (trim(arg)=="-alpha") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') alpha
            ci = ci + 2
        else if (trim(arg)=="-beta") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') beta
            ci = ci + 2
        else if (trim(arg)=="-dt") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') dt
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
    
    ! setting parameters for grids and calculations
    m = 1.0d0                    ! mass of the electron
    ndim = n**2                  ! size of 2e- direct product space
    dx = 2.0d0*x0 / real(n-1,8)    ! grid-spacing
    multi_well_switch = .false.  ! default for single well
    te_swtich = .true.           ! switch to test one-electron hamiltonian 
    
    print*, ""
    print*, "    Two-electron DVR    "
    print*, "------------------------"
    print('(a,f16.8)'), "x0     = ", x0
    print('(a,f16.8)'), "m      = ", m
    print('(a,f16.8)'), "Z      = ", z
    print('(a,i16)'),   "n      = ", n
    print('(a,i16)'),   "ndim   = ", ndim
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
    
    nsparse = 2*ndim**2 - 1
    ! Allocating arrays
    allocate(h_array(nsparse))
    allocate(h_row(ndim+1))
    allocate(h_col(nsparse))
    allocate(psi0(ndim))

    h_array = 0.0d0
    print*, "Generating Hamiltonian Matrix......"
    print*, "" 
    call sparse_te_sw_hamiltonian(h_array, h_row, h_col)
    print*, "Completed generating Hamiltonian Matrix!"
    print*, "" 
    k = size(h_array)
    
    print*, "Using Imaginary time propagation (ITP) to get ground state......"

    call gen_trial_state(psi0)
    allocate(psi(ndim))
    call itp_sparse(h_array, h_row, h_col, psi0, dt, print_nstep, etol, Ei, tstep)
    comment = "! imag_tprop final wavefunction"
    psi = cmplx(psi0, 0.0d0, kind=dp)
    deallocate(psi0)
    open(100, file='psi_itp.wfn')
        call write_wfn(100, comment, psi)
    close(100)
    print*, ""
    print*, "Completed Imaginary time propagation (ITP)."
    print*, ""
    print*, "ITP Summary:"
    print*, "------------"
    print*, "dt      = ", dt
    print*, "nstep   = ", tstep
    print*, "E_conv  = ", etol
    print*, "Ei      = ", Ei
    print*, ""
end program main
