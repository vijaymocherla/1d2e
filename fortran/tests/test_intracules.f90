program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use intracules
    use helpers, only: read_wfn, write_wfn
    implicit none
    
    complex(dp), allocatable, dimension(:) :: wfn
    real(dp), allocatable, dimension(:) :: wigner_intracule
    real(dp) :: x0
    integer :: i, m
    character (len=128) :: wfn_file
    character (len=128) :: out_file  
    character (len=128) :: comment
    integer :: ci
    character (len=128) :: arg 
    
    wfn_file='test.wfn'
    out_file = 'wigner_test.txt'
    ci = 1
    do 
        call get_command_argument(ci, arg)
        if (trim(arg)=="-i") then
            call get_command_argument(ci+1,arg)
            read(arg, '(a128)') wfn_file
            ci = ci + 2
        else if (trim(arg)=="-o") then
            call get_command_argument(ci+1,arg)
            read(arg, '(a128)') out_file
            ci = ci + 2
        else
            exit
        end if
    end do

    open(100,file=trim(wfn_file))
    call read_wfn(100, comment, wfn)
    close(100)
    x0 = 10.0d0
    call comp_wigner_intracule(wfn, x0, wigner_intracule)
    open(300, file=trim(out_file))
    m = size(wigner_intracule)
    write(300,*) m
    write(300,*) trim(comment)
    do i=1,m
            write(300,*) wigner_intracule(i)
        end do
    close(300)
    
end program main