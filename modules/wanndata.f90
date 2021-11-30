!
!   wanndata.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!

MODULE wanndata
  !
  use constants
  !
  IMPLICIT NONE

  INTEGER norb
  INTEGER nrpt
  INTEGER r000 ! This is the on-site

  COMPLEX(DP), ALLOCATABLE :: ham(:,:,:)

  REAL(DP), ALLOCATABLE :: weight(:)

  REAL(DP), ALLOCATABLE :: rvec(:,:)
  !
CONTAINS

SUBROUTINE ham_shift_ef(mu)
  !
  USE constants
  implicit none
  !
  real(dp) :: mu
  !
  integer ii
  !
  do ii=1, norb
    ham(ii, ii, r000) = ham(ii, ii, r000)-mu
  enddo
  !
END SUBROUTINE

SUBROUTINE read_ham(seed)
!
  USE para,      only: para_sync_int, para_sync_real, para_sync_cmplx, inode
  USE constants
  !
  IMPLICIT NONE
  !
  CHARACTER(len=80) seed
  INTEGER irpt, iorb, jorb
  integer, dimension(5) :: tt
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  if (inode.eq.0) then
    write(stdout, *) " # Reading file "//trim(seed)//"_hr.dat"
    !
    open(unit=fin, file=trim(seed)//"_hr.dat")
    !
    read(fin, *)
    read(fin, *) tt(1) ! norb
    read(fin, *) tt(2) ! nrpt
    !
    write(stdout, *) " #  Dimensions:"
  endif
  !
  CALL para_sync_int(tt, 2)
  norb=tt(1)
  nrpt=tt(2)
  !
  allocate(ham(1:norb, 1:norb, 1:nrpt))
  allocate(weight(1:nrpt))
  allocate(rvec(1:3, 1:nrpt))
  !
  if (inode.eq.0) then
    write(stdout, *) "    # of orbitals:", norb
    write(stdout, *) "    # of real-space grid:", nrpt
    allocate(wt(1:nrpt))
    read(fin, '(15I5)') (wt(irpt),irpt=1,nrpt)
    weight(:)=wt(:)
    deallocate(wt)
    !
    do irpt=1, nrpt
      do iorb=1, norb
        do jorb=1, norb
          read(fin, *) tt, a, b
          if ((jorb.eq.1).and.(iorb.eq.1)) then
            rvec(:, irpt)=tt(1:3)
          endif
          ham(jorb, iorb, irpt)=CMPLX(a, b, KIND=dp)
        enddo
      enddo
    enddo
    !
    close(unit=fin)
    write(stdout, *) " # Done."
  endif
  !
  CALL para_sync_cmplx(ham, norb*norb*nrpt)
  CALL para_sync_real(weight, nrpt)
  CALL para_sync_real(rvec, 3*nrpt)
  !
  do irpt=1, nrpt
    if (sum(rvec(:, irpt)**2)<eps4) r000=irpt
  enddo
  !
END SUBROUTINE

SUBROUTINE finalize_wann()
  !
  IMPLICIT NONE
  !
  if (allocated(ham)) deallocate(ham)
  if (allocated(weight)) deallocate(weight)
  if (allocated(rvec)) deallocate(rvec)
  !
END SUBROUTINE

END MODULE

