MODULE banddata
  !
  use constants
  !
  implicit none
  !
  real(dp) ef
  real(dp), allocatable :: eig(:, :)
  complex(dp), allocatable :: egv(:, :, :)
  integer nbnd
  integer nkx, nky, nkz
  integer nkpt  ! = nkx*nky*nkz
  !
CONTAINS

SUBROUTINE finalize_band()
  !
  implicit none
  !
  if(allocated(eig)) deallocate(eig)
  if(allocated(egv)) deallocate(egv)
  !
END SUBROUTINE

END MODULE

