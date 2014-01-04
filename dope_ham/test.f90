PROGRAM test
  !
  use constants
  use orbitals
  use dopedata
  !
  IMPLICIT NONE
  !
  real(dp) tvec(1:3)
  integer ii
  norb=40
  tvec(1)=0.50
  tvec(2)=0.00
  tvec(3)=0.00
  CALL read_orbitals
  CALL initialize_dopedata
  CALL translate_crystal(tvec)
  do ii=1, norb
    write(*,*) new_to_old(ii), new_to_old_dR(ii, :)
  enddo
  CALL finalize_orbitals
  !
END PROGRAM
