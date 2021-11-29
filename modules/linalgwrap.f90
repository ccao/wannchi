MODULE linalgwrap
  !
  use constants, only : dp
  !
  implicit none
  !
  interface invmat
    module procedure dinvmat, zinvmat
  end interface
  !
  contains
  !
  subroutine dinvmat(xmat, ndim)
    !
    use constants, only : dp
    !
    implicit none
    !
    integer :: ndim
    real(dp), dimension(ndim, ndim) :: xmat
    !
    real(dp), dimension(ndim) :: work
    integer, dimension(ndim) :: ipiv
    integer :: info
    !
    call dgetrf(ndim, ndim, xmat, ndim, ipiv, info)
    call dgetri(ndim, xmat, ndim, ipiv, work, ndim, info)
    !
  end subroutine
  !
  subroutine zinvmat(xmat, ndim)
    !
    use constants, only : dp
    !
    integer :: ndim
    complex(dp), dimension(ndim, ndim) :: xmat
    !
    complex(dp), dimension(ndim) :: work
    integer, dimension(ndim) :: ipiv
    integer :: info
    !
    call zgetrf(ndim, ndim, xmat, ndim, ipiv, info)
    call zgetri(ndim, xmat, ndim, ipiv, work, ndim, info)
    !
  end subroutine
  !
END MODULE

