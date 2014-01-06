MODULE banddata
  !
  use constants
  !
  implicit none
  !
  real(dp) ef
  real(dp), allocatable :: occ(:, :)
  real(dp), allocatable :: eig(:, :)
  complex(dp), allocatable :: egv(:, :, :)
  real(dp), allocatable :: kvec(:, :)
  integer nbnd
  integer nkx, nky, nkz
  integer nkpt  ! = nkx*nky*nkz = dimension of kvec
  !
CONTAINS

SUBROUTINE init_band
  !
  implicit none
  !
  integer ikx, iky, ikz, ik
  !
  allocate(occ(1:nbnd, 1:nkpt))
  allocate(eig(1:nbnd, 1:nkpt))
  allocate(egv(1:nbnd, 1:nbnd, 1:nkpt))
  allocate(kvec(1:3, 1:nkpt))
  !
  occ(:,:)=0.d0
  eig(:,:)=0.d0
  egv(:, :, :)=cmplx_0
  !
  do ikx=1, nkx
    do iky=1, nky
      do ikz=1, nkz
        ik=(ikz-1)*nkx*nky+(iky-1)*nkx+ikx
        kvec(1, ik)=(ikx-1.d0)/nkx
        kvec(2, ik)=(iky-1.d0)/nky
        kvec(3, ik)=(ikz-1.d0)/nkz
      enddo
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE finalize_band()
  !
  implicit none
  !
  if(allocated(occ)) deallocate(occ)
  if(allocated(eig)) deallocate(eig)
  if(allocated(egv)) deallocate(egv)
  if(allocated(kvec)) deallocate(kvec)
  !
END SUBROUTINE

END MODULE

