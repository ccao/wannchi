MODULE dopedata
  !
  use constants
  !
  IMPLICIT NONE

  INTEGER ndopant, ref_dopant
 
  INTEGER, ALLOCATABLE :: new_to_old(:)

  INTEGER, ALLOCATABLE :: new_to_old_dR(:, :)

  INTEGER, ALLOCATABLE :: dopant_list(:)

  COMPLEX(DP), ALLOCATABLE :: doped_ham(:,:,:)  !

  COMPLEX(DP), ALLOCATABLE :: dham(:,:,:)
  !
CONTAINS
  !
  SUBROUTINE initialize_dopedata
    !
    use constants
    use orbitals, only: nat
    use wanndata, only: norb, nrpt
    !
    IMPLICIT NONE
    !
    allocate(new_to_old(1:nat))
    allocate(new_to_old_dR(1:nat, 1:3))
    allocate(dopant_list(1:ndopant))
    allocate(doped_ham(1:norb, 1:norb, 1:nrpt))
    allocate(dham(1:norb, 1:norb, 1:nrpt))
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_dopedata
    !
    IMPLICIT NONE
    !
    if(allocated(new_to_old)) deallocate(new_to_old)
    if(allocated(new_to_old_dR)) deallocate(new_to_old_dR)
    if(allocated(dopant_list)) deallocate(dopant_list)
    if(allocated(doped_ham)) deallocate(doped_ham)
    if(allocated(dham)) deallocate(dham)
    !
  END SUBROUTINE
  !
END MODULE

