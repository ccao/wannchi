PROGRAM wannband
  !
  use para
  use constants
  use wanndata, only: finalize_wann
  use banddata, only: eig, nbnd, finalize_band, sigma
  use input, only: code, nqseg, nqbnd, get_qvec
  !
  implicit none
  !
  real(dp) qv(1:3)
  integer iq
  character(len=80) seed
  !
  code=2
  seed="wannier90"
  !
  CALL init_para
  !
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL read_input
  !
  do iq=1, nqseg*nqbnd+1
    CALL get_qvec(qv, iq)
  enddo
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM
