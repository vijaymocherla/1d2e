! 
!   Real time propagation of 2e- wavefunction
! 
program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use iso_c_binding
    use lapack_wrappers, only:eigsh
    use blas_wrappers
    use real_tprop

    implicit none
    complex(dp), allocatable, target :: h_array(:)
    real(dp), pointer :: h_array_real(:)=> null()
    integer,  allocatable :: h_row(:)
    integer,  allocatable :: h_col(:)
    complex(dp), allocatable :: psi0(:)
    integer :: ndim              ! size of 2e- matrix
    integer :: nsparse           ! no. of non-zero elements per row
    integer  :: i, tstep         ! iteration variables
    integer  :: print_nstep      ! no. of time-steps after which to print
    real(dp) :: Ei               ! energy at ti
    real(dp) :: dx1, dx2         ! shift in wavefxn position
    real(dp) :: dt               !  
    real(dp) :: ti, tf           !
    integer  :: ci               ! 
    character (len=32) :: arg    !
    
    ci = 1
    ! Reading input parameters
    n = 64                       ! size of 1e- grids
    x0 = 10.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    ti = 0.0d0                   ! initial time 
    tf = 1.000d0                 ! final time
    dt  = 0.001d0                ! time step of propagation
    print_nstep = 100
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
        else if (trim(arg)=="-ti") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') ti
            ci = ci + 2    
        else if (trim(arg)=="-tf") then
            call get_command_argument(ci+1,arg)
            read(arg, '(f32.16)') tf
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
    nsparse = 2*n - 1
    allocate(h_array(ndim*nsparse))
    allocate(h_col(ndim*nsparse))
    allocate(h_row(ndim+1))
    allocate(psi0(ndim))

    h_array = cmplx(0.0d0,0.0d0, kind=dp)
    print*, "Generating Hamiltonian Matrix......"
    print*, "" 
    ! Some programming stunts using pointers to avoid creating a new array.
    ! Assuming h_array is contigous in memory, assigning a c_pointer of size=ndim
    ! to a fortran_pointer (h_array_real). This effectively means assigning 
    ! first ndim addresses of complex array (h_array) to the fortran_pointer.
    call c_f_pointer(c_loc(h_array), h_array_real, shape=[size(h_array)])
    call sparse_te_sw_hamiltonian(h_array_real, h_row, h_col)
    ! Assigning h_array_real addresses real part of complex array and setting 
    ! imag. part to zero.    
    h_array%re = h_array_real
    h_array%im = 0.0d0

    print*, "Completed generating Hamiltonian Matrix!"
    print*, "" 
    
    print*, "Propagating the wavefunction......"

    dx1 = 0.0d0
    dx2 = 0.0d0
    call gen_intial_wfn(psi0, dx1, dx2)
    call rtp_sparse(h_array, h_row, h_col, psi0, ti, tf, dt, print_nstep, Ei, tstep)
    open(100, file='final.wfn')
    do i=1,ndim
        write(100,*) psi0(i)
    end do
    close(100)
    print*, ""
    print*, "   Summary of Real time propagation   "
    print*, "--------------------------------------"
    print*, "n       = ", n
    print*, "dt      = ", dt
    print*, "ti      = ", ti
    print*, "tf      = ", tf
    print*, "Ei      = ", Ei
    print*, ""
end program main
