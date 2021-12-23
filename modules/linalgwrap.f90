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
  subroutine sparsemulmat(zmat, xmat_cp, ymat, idxcp, ndim, nidxcp, alpha, beta)
    !
    use constants, only : dp
    !
    implicit none
    !
    integer :: ndim
    ! Dimension of full matrices:
    !  xmat(ndim, ndim); ymat(ndim, ndim)
    integer :: nidxcp
    ! Dimension of nonzero elements in ymat
    integer, dimension(2, nidxcp) :: idxcp
    ! xmat(idxcp(1, ii), idxcp(2, ii))=xmat_cp(ii)
    complex(dp), dimension(ndim, ndim) :: ymat, zmat
    ! input & output matrices
    real(dp), dimension(nidxcp) :: xmat_cp
    real(dp) :: alpha, beta
    ! compact form of ymat
    integer ii, jj, i1, i2
    !
    zmat(:, :)=beta*zmat(:, :)
    !
    do ii=1, nidxcp
      i1=idxcp(1, ii)
      i2=idxcp(2, ii)
      do jj=1, ndim
        zmat(i1, jj)=zmat(i1, jj)+alpha*xmat_cp(ii)*ymat(i2, jj)
      enddo
    enddo
    !
  end subroutine

  subroutine matmulsparse(zmat, xmat, ymat_cp, idxcp, ndim, nidxcp, alpha, beta)
    !
    use constants, only : dp, cmplx_0
    !
    implicit none
    !
    integer :: ndim
    ! Dimension of full matrices:
    !  xmat(ndim, ndim); ymat(ndim, ndim)
    integer :: nidxcp
    ! Dimension of nonzero elements in ymat
    integer, dimension(2, nidxcp) :: idxcp
    ! ymat(idxcp(1, ii), idxcp(2, ii))=ymat_cp(ii)
    complex(dp), dimension(ndim, ndim) :: xmat, zmat
    ! input & output matrices
    real(dp), dimension(nidxcp) :: ymat_cp
    real(dp) :: alpha, beta
    ! compact form of ymat
    !
    integer ii, jj, j1, j2
    !
    zmat(:, :)=beta*zmat(:, :)
    do ii=1, ndim
      do jj=1, nidxcp
        j1=idxcp(1, jj)
        j2=idxcp(2, jj)
        zmat(ii, j2)=zmat(ii, j2)+alpha*xmat(ii, j1)*ymat_cp(jj)
      enddo
    enddo
    !
  end subroutine

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

