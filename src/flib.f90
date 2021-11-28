INTEGER FUNCTION kpt_index(kv)
  !
  use banddata, only : nkx, nky, nkz
  use constants, only : dp, eps6
  !
  implicit none
  !
  real(dp) kv(1:3)
  integer ik, ikx, iky, ikz
  !
  ikx=nint((kv(1)-nint(kv(1)))*nkx)
  iky=nint((kv(2)-nint(kv(2)))*nky)
  ikz=nint((kv(3)-nint(kv(3)))*nkz)
  !
  kpt_index=ikz*nkx*nky+iky*nkx+ikx+1
  return
  !
END FUNCTION
