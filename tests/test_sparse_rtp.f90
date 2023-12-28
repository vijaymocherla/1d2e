! 
!   Real time propagation of 2e- wavefunction
! 
program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use two_electron_dvr
    use iso_c_binding
    use blas_wrappers
    use real_tprop
    use helpers, only: write_wfn

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
    real(dp) :: p01              !
    real(dp) :: p02              ! shift in wavefxn position
    real(dp) :: zeta             ! effective nuclear charge in efn ansatz
    real(dp) :: dt               !  
    real(dp) :: ti, tf           !
    integer  :: ci, rc           ! 
    character (len=32) :: arg    !
    character (len=32) :: input_file
    character (len=32) :: output_file
    character (len=128):: comment

    namelist /inp_rtp/ n, x0, z, alpha, beta, &
                          multi_well_switch, te_swtich, & 
                          p01, p02, zeta, &
                          ti, tf, dt, print_nstep 
        
    ! Default input parameters
    n = 32                       ! size of 1e- grids
    x0 = 15.0d0                  ! extent of 1d box
    alpha = 1.00d0               ! 1e- soft coulomb parameter
    beta  = 1.00d0               ! e- correlation parameter
    ti = 0.0d0                   ! initial time 
    tf = 1.000d0                 ! final time
    dt  = 0.001d0                ! time step of propagation
    print_nstep = 100
    p01 = 0.0d0
    p02 = 0.0d0
    zeta = 1.0d0 
    ! setting parameters for grids and calculations
    z = 2.0d0
    m = 1.0d0                    ! mass of the electron
    multi_well_switch = .false.  ! default for single well
    te_swtich = .true.           ! switch to test one-electron hamiltonian 
    ci = 1
    output_file = 'real_tprop.out' 
    ! Checking for input file 
    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-i") then
            call get_command_argument(ci+1,arg)
            read(arg, '(a32)') input_file
            ci = ci + 2
            ! reading input parameters
            open(100, file=trim(input_file))
            read(100, nml=inp_rtp, iostat=rc)
            if (rc /= 0) then
                print*, '!! I/O ERROR: Can not read input file. IOSTAT=',rc
                stop
            end if
            close(100)
        else if (trim(arg)=="-o") then
            call get_command_argument(ci+1,arg)
            read(arg, '(a32)') output_file
            ci = ci + 2            
        else
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
        end if        
    end do
    
        
    ndim = n**2                  ! size of 2e- direct product space
    dx = 2.0d0*x0 / real(n-1,8)  ! grid-spacing
        
    open(200, file=trim(output_file))
    write(200,*) ""
    write(200,*) "    Two-electron DVR    "
    write(200,*) "------------------------"
    write(200,'(a,f16.8)') "x0     = ", x0
    write(200,'(a,f16.8)') "m      = ", m
    write(200,'(a,f16.8)') "Z      = ", z
    write(200,'(a,i16)')   "n      = ", n
    write(200,'(a,i16)')   "ndim   = ", ndim
    write(200,'(a,f16.8)') "dx     = ", dx
    write(200,*) "" 
    write(200,'(a)') "Soft-Coulomb Parameters:"
    write(200,'(a,f16.8)') "alpha  = ", alpha
    write(200,'(a,f16.8)') "beta   = ", beta
    write(200,*) "" 
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
    write(200,*) "Generating Hamiltonian Matrix......"
    write(200,*) "" 
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
    write(200,*) "Completed generating Hamiltonian Matrix!"
    write(200,*) "" 
    
    write(200,*) "    Initial wavefxn ansatz    "
    write(200,*) "------------------------------"
    write(200,'(a,f16.8)') "zeta    = ", zeta
    write(200,'(a,f16.8)') "p01     = ", p01
    write(200,'(a,f16.8)') "p02     = ", p02
    call gen_intial_wfn(psi0, p01, p02, zeta)
    write(200,*) ""

    write(200,*) ""
    write(200,*) "        Real time propagation         "
    write(200,*) "--------------------------------------"
    write(200,*) "n       = ", n
    write(200,*) "dt      = ", dt
    write(200,*) "ti      = ", ti
    write(200,*) "tf      = ", tf
    write(200,*) ""
    write(200,*) "Propagating the wavefunction......"
    call rtp_sparse(h_array, h_row, h_col, psi0, ti, tf, dt, print_nstep, Ei, tstep)
    write(200,*) ""
    write(200,*) "Completed Real-time propagation!"
    write(200,*) ""
    ! writing final wavefxn to txt file
    open(100, file='psi_tf_rtp.wfn')
    write(comment,'(a,f16.8)') "! final wavefunction from sparse real-time propagation, at tf =",tf
    call write_wfn(100, comment, psi0)
    close(100)
end program main
