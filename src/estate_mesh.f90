MODULE estate_mesh
  !
  use constants
  !
  implicit none
  !
  real(dp), allocatable :: eig(:, :)
  complex(dp), allocatable :: egv(:, :, :)
  integer nkx, nky, nkz
  integer nkmesh
  !
END MODULE
