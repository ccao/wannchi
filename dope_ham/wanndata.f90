!
!   wanndata.f90
!   
!
!   Created by Chao Cao on 1/22/13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!

MODULE wanndata
  use constants
  IMPLICIT NONE

  LOGICAL first_ham

  INTEGER norb
  INTEGER nrpt

  COMPLEX(DP), ALLOCATABLE :: ham(:,:,:,:)  ! ham(0, :, :, :)  & ham(1, :, :, :)

  INTEGER, ALLOCATABLE :: weight(:)

  INTEGER, ALLOCATABLE :: rvec(:,:)

  REAL(DP) ef
CONTAINS

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

