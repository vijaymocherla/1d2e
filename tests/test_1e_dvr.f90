program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use lapack_wrappers, only:eigsh
    use imag_tprop
    implicit none
    real(dp), allocatable :: hamiltonian(:,:), evecs(:,:)  
    real(dp), allocatable :: evals(:)
    real(dp), allocatable :: psi0(:), psi(:)
    integer :: i,j, tstep, ndim
    real(dp) :: dt
    real(dp) :: etol, Ei
    

    ! setting parameters for grids
    m = 1.0d0                   ! mass of the electron
    n = 100                     ! size of 1e- grids
    ndim = n                    ! size of 2e- direct product space
    x0 = 10.0                   ! extent of 1d box
    dx = 2.0d0*x0 / real(n-1,8) ! grid-spacing
    alpha = 1.00d0              ! 1e- soft coulomb parameter
    multi_well_switch = .false. ! default for single well
    te_swtich = .false.         ! switch to test one-electron hamiltonian 
    etol = 0.000000001          ! energy absolute tolerance 

    print*, ""
    print*, "     1e- Single-well DVR     "
    print*, "-----------------------------"
    print('(a,f16.8)'), "x0     = ", x0
    print('(a,i16)'),   "n      = ", n
    print('(a,i16)'),   "ndim   = ", ndim
    print('(a,f16.8)'), "dx     = ", dx
    print('(a,f16.8)'), "alpha  = ", alpha
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
    allocate(psi0(ndim), psi(ndim))

    hamiltonian = 0.0d0
    print*, "Generating Hamiltonian Matrix......"
    print*, "" 
    do i=1,ndim
        hamiltonian(i,i) = oe_scsw_hamiltonian(i,i)
        do j=i+1,ndim
            hamiltonian(i,j) = oe_scsw_hamiltonian(i,j)
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
    open(100, file='evec0.txt')
    do i=1,ndim
        write(100,*) evecs(i,1)
    end do
    close(100)

    
    print*, "Using Imaginary time propagation to get ground state......"
    dt  = 0.0001d0

    call gen_trial_state(psi0)    

    call itp_on_the_fly(psi0, dt, 100, etol, Ei, tstep)
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
