SUBROUTINE map_to_kpt(kv, ik, nkx, nky, nkz)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) kv(1:3)
  integer ik, nkx, nky, nkz
  !
  integer ikx, iky, ikz
  !
  ikx=mod(ik-1, nkx)
  iky=mod( (ik-1-ikx)/nkx, nky )
  ikz=(ik-1-iky*nkx-ikx)/(nkx*nky)
  kv(1)=ikx*(1.d0/nkx)
  kv(2)=iky*(1.d0/nky)
  kv(3)=ikz*(1.d0/nkz)
  !
END SUBROUTINE

SUBROUTINE map_to_idx(ik, kv, nkx, nky, nkz)
  !
  use constants, only: dp
  !
  implicit none
  !
  integer ik, nkx, nky, nkz
  !
  real(dp) kv(1:3)
  !
  integer ikx, iky, ikz, iik
  !
  ikx=kv(1)*nkx
  iky=kv(2)*nky
  ikz=kv(3)*nkz
  !
  ikx=mod(ikx, nkx)
  iky=mod(iky, nky)
  ikz=mod(ikz, nkz)
  !
  if (ikx<0) ikx=ikx+nkx
  if (iky<0) iky=iky+nky
  if (ikz<0) ikz=ikz+nkz
  !
  ik=ikz*nky*nkx+iky*nkx+ikx+1
  !
END SUBROUTINE
