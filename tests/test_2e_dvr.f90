program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use lapack_wrappers, only:eigsh
    use imag_tprop
    use helpers, only: write_wfn

    implicit none
    real(dp), allocatable :: hamiltonian(:,:), evecs(:,:)  
    real(dp), allocatable :: evals(:)
    complex(dp), allocatable :: psi(:)
    real(dp), allocatable :: psi0(:)
    integer :: i,j, tstep, print_nstep, ndim
    integer :: print_nevecs
    real(dp) :: dt
    real(dp) :: etol, Ei
    integer :: ci
    character (len=32) :: arg
    character (len=128) :: comment
    
    ci = 1
    ! Reading input parameters
    n = 8                        ! size of 1e- grids
    x0 = 10.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    etol = 0.000000001           ! energy absolute tolerance 
    dt  = 0.001d0
    print_nstep = 100
    print_nevecs = 20            ! print first n eigen-vectors
    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-n") then
            call get_command_argument(ci+1,arg)
            read(arg, '(i32)') n
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
        else if (trim(arg)=="-print_nstep") then
            call get_command_argument(ci+1,arg)
            read(arg, '(i32)') print_nstep
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
    

    ! Allocating arrays
    allocate(hamiltonian(ndim, ndim))
    allocate(evals(ndim))
    allocate(evecs(ndim, ndim))
    allocate(psi0(ndim))
    allocate(psi(ndim))
    
    hamiltonian = 0.0d0
    print*, "Generating Hamiltonian Matrix......"
    print*, "" 
    do i=1,ndim
        hamiltonian(i,i) = te_sw_hamiltonian(i,i)
        do j=i+1,ndim
            hamiltonian(i,j) = te_sw_hamiltonian(i,j)
            hamiltonian(j,i) = hamiltonian(i,j)
        end do
    end do
   
    print*, "Diagonalising the Hamiltonian......"
    print*, "" 
    ! diagonalising hamiltonian
    call eigsh(hamiltonian, evals, evecs)
    print*, "Completed Diagonalisation!"
    print*, "Ground state Energy: ", evals(1)
    print*,""
    open(100, file='evals.txt')
    do i=1,ndim
        write(100,*) evals(i)
    end do
    close(100)
    open(100, file='evecs.txt')
    do i=1,ndim
        write(100,*) evecs(i,1:print_nevecs)
    end do

    
    print*, "Using Imaginary time propagation (ITP) to get ground state......"

    call gen_trial_state(psi0)   
    call itp_on_the_fly(psi0, dt, print_nstep, etol, Ei, tstep)
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
