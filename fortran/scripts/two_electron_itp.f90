program main
    use two_electron_dvr
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use lapack_wrappers, only:eigsh
    use blas_wrappers
    implicit none
    real(dp), allocatable :: psi0(:)
    integer :: i, tstep, print_nstep
    real(dp) :: dt
    real(dp) :: etol, Ei
    integer :: ci
    character (len=32) :: arg
    character (len=32) :: input_file
    character (len=32) :: output_file
    namelist /itp_1d2e/ n, x0, alpha, beta, dt, print_nstep, etol
    ci = 1
    ! Default parameters
    n = 32                       ! size of 1e- grids
    x0 = 10.0d0                  ! extent of 1d box
    z  = 2.0d0                   ! atomic number
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    
    etol = 0.000000001           ! energy absolute tolerance 
    dt  = 0.001d0
    print_nstep = 100

    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-i") then
            call get_command_argument(ci+1,arg)
            read(arg, *) input_file
            ci = ci + 2
        else if (trim(arg)=="-o") then
            call get_command_argument(ci+1,arg)
            read(arg, *) output_file
            ci = ci + 2
        else
            exit
        end if        
    end do
    

    ! setting parameters for grids and calculations
    z = 2.0d0
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
    
    
    print*, "Using Imaginary time propagation (ITP) to get ground state......"

    call gen_trial_state(psi0)   
    call imag_tprop(psi0, dt, print_nstep, etol, Ei, tstep)
    open(100, file='psi_itp.out')
    do i=1,ndim
        write(100,*) psi0(i)
    end do
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