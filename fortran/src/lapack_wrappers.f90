module lapack_wrappers
    use, intrinsic ::iso_fortran_env, only:dp=>real64
    implicit none
    private

    public eigsh 

    contains

    ! wrapper for DSYEV
    subroutine eigsh(A, vals, vecs)
        real(dp), dimension(:, :), intent(in) :: A
        real(dp), dimension(:), intent(out) :: vals
        real(dp), dimension(:, :), intent(out) :: vecs
        integer(dp) :: info, ndim
        real(dp), allocatable :: work(:)
        character, parameter :: jobz="V", uplo="U"
        ndim = size(A, dim=1) 
        allocate(work(3*ndim))
        vecs = A
        call dsyev(jobz, uplo, ndim, vecs, ndim, vals, work, 3*ndim, info)
    end subroutine eigsh

end module lapack_wrappers