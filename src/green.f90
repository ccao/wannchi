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

