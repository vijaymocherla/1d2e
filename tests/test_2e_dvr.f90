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
    integer :: i,j, ndim
    integer :: print_nevecs
    integer :: ci
    character (len=32) :: arg
    character (len=128) :: comment
    
    ci = 1
    ! Reading input parameters
    n = 8                        ! size of 1e- grids
    x0 = 15.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    z = 2.0d0                    ! atomic number   
    print_nevecs = 20            ! print first n eigen-vectors

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
    

    ! Allocating arrays
    allocate(hamiltonian(ndim, ndim))
    allocate(evals(ndim))
    allocate(evecs(ndim, ndim))
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
    deallocate(hamiltonian)
    print*, "Completed Diagonalisation!"
    print*,""
    print*, "Saving all eigen-values to 'evals.txt'"
    print*,""
    open(100, file='evals.txt')
    do i=1,ndim
        write(100,*) evals(i)
    end do
    close(100)
    print('(a,i4,a)'), "Saving first ",print_nevecs," eigen-vectors to 'evecs.txt'"
    print*,""
    open(100, file='evecs.txt')
    do i=1,ndim
        write(100,*) evecs(i,1:print_nevecs)
    end do
    close(100)
    print*, "Ground state Energy (in Eh) : ", evals(1)
    print*,""
    ! converting (real) psi0 to (complex) psi
    psi = cmplx(evecs(:,1), 0.0d0, kind=dp)
    comment = "! exact ground state wavefunction"
    open(200, file='psi_exact.wfn')
    call write_wfn(200, comment, psi)
    close(200)
    deallocate(evals, evecs, psi)
    
end program main
