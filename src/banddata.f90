MODULE banddata
  !
  use constants
  !
  implicit none
  !
  real(dp) ef
  real(dp), allocatable :: eig(:, :)
  complex(dp), allocatable :: egv(:, :, :)
  complex(dp), allocatable :: xi(:, :, :, :)
  real(dp), allocatable :: kmesh(:, :)
  real(dp), allocatable :: qvec(:, :)
  integer nbnd
  integer nkx, nky, nkz
  integer nkpt  ! = nkx*nky*nkz = dimension of kmesh
  integer nqpt  ! = dimension of qvec
  real(dp) sigma
  !
CONTAINS

SUBROUTINE finalize_band()
  !
  implicit none
  !
  if(allocated(eig)) deallocate(eig)
  if(allocated(egv)) deallocate(egv)
  if(allocated(kmesh)) deallocate(kmesh)
  if(allocated(qvec)) deallocate(qvec)
  !
END SUBROUTINE

END MODULE

