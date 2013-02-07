PROGRAM wannchi
  !
  use constants
  use wanndata, only: norb, finalize_wann
  use banddata, only: nbnd, nkx, nky, nkz, finalize_band
  !
  implicit none
  !
  integer nqpt
  integer iq
  real(dp) qvec(1:3)
  complex(dp) chi
  character(len=80) seed
  !
  seed="wannier90"
  !
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL read_input
  !
  CALL interpolate_bands
  !
  nqpt=nkx*nky
  !
  do iq=1, nqpt
    !
    CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
    CALL compute_chi_diag(chi, qvec)
    !
    write(*,'(3F12.8,2F22.12)') qvec, chi
    !
    if (abs(qvec(1)).le.eps) write(*,*) ' '
  enddo
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
END PROGRAM
