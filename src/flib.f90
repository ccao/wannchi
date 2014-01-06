INTEGER FUNCTION kpt_index(kv)
  !
  use banddata, only : kvec, nkpt
  use constants, only : dp, eps6
  !
  implicit none
  !
  real(dp) kv(1:3)
  integer ik
  logical found
  real(dp) dk(1:3)
  !
  do ik=1, nkpt
    dk(:)=(kv(:)-kvec(:, ik))-nint(kv(:)-kvec(:, ik))
    if ( abs(dk(1))<eps6.and.abs(dk(2))<eps6.and.abs(dk(3))<eps6 ) then
      kpt_index=ik
      return
    endif
  enddo
  !
  kpt_index=-1
  return
  !
END FUNCTION
