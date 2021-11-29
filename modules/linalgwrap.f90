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
  interface eigen
    module procedure heigen, geigen
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
  subroutine heigen(eig, xmat, ndim)
    ! Symmetric (hermitian) eigen value problem
    !
    use constants, only : dp
    !
    integer :: ndim
    complex(dp), dimension(ndim, ndim) :: xmat
    real(dp), dimension(ndim) :: eig
    !
    integer info
    complex(dp), dimension(2*ndim) :: work
    complex(dp), dimension(3*ndim) :: rwork
    !
    ! call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    call zheev('V', 'U', ndim, xmat, ndim, eig, work, 2*ndim, rwork, info)
    !
  end subroutine
  !
  subroutine geigen(eig, xmat, ndim)
    ! Generic eigen value problem
    ! Solves right eigen states
    !
    use constants, only : dp
    !
    integer :: ndim
    complex(dp), dimension(ndim, ndim) :: xmat
    complex(dp), dimension(ndim) :: eig
    !
    integer info
    complex(dp), dimension(2*ndim) :: work
    complex(dp), dimension(2*ndim) :: rwork
    complex(dp), dimension(1, 1)   :: vl
    complex(dp), dimension(ndim, ndim) :: vr
    !
    !call zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    call zgeev('N', 'V', ndim, xmat, ndim, eig, vl, 1, vr, ndim, work, 2*ndim, rwork, info)
    xmat(:, :)=vr(:, :)
    !
  end subroutine
  !
END MODULE

