SUBROUTINE calc_g0(gf, hk, w, ndim, inv)
  !
  use constants,       only: dp
  use linalgwrap,      only: invmat
  !
  implicit none
  !
  integer ndim
  complex(dp), dimension(ndim, ndim) :: gf, hk
  complex(dp) :: w
  logical :: inv ! if true, gf contains G^-1 instead of G
  !
  integer ii
  !
  gf(:,:)=-hk(:,:)
  do ii=1, ndim
    gf(ii, ii)=w+gf(ii, ii)
  enddo
  !
  if (.not.inv) then
    call invmat(gf, ndim)
  endif
  !
END SUBROUTINE

SUBROUTINE calc_corrFF_gf(gf, hff, sigmat, Ecc_diag, Vfc, w, ndimf, ndimc, inv)
  !
  use constants,       only: dp, cmplx_0, cmplx_1
  use linalgwrap,      only: invmat
  !
  implicit none
  !
  integer ndimf, ndimc
  real(dp), dimension(ndimc) :: Ecc_diag
  complex(dp), dimension(ndimf, ndimf) :: gf, hff, sigmat
  complex(dp), dimension(ndimf, ndimc) :: Vfc
  complex(dp)   :: w
  logical       :: inv
  !
  integer ii, jj
  complex(dp), dimension(ndimc) :: gcc_diag
  !
  !gf(:, :)=cmplx_0
  gcc_diag(:)=cmplx_1/(w-Ecc_diag(:))
  !
  ! Hybridization = Vfc * gcc * Vfc^{\dagger}
  do ii=1, ndimf
    do jj=1, ndimf
      gf(ii, jj)=-sum(Vfc(ii, :)*gcc_diag(:)*conjg(Vfc(jj, :)))
    enddo
  enddo
  !
  !write(*, *) " -Delta:"
  !call show_matrix(gf, ndimf, stdout)
  !
  gf=gf-hff-sigmat
  !
  do ii=1, ndimf
    gf(ii, ii)=w+gf(ii, ii)
  enddo
  !
  if (.not.inv) then
    call invmat(gf, ndimf)
  endif
  !
END SUBROUTINE
