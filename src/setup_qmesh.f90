SUBROUTINE setup_kmesh
  !
  use para
  use input    : only mode, nqseg, nqbnd, get_qvec
  use banddata : only qvec
  !
  implicit none
  !
  integer nq, iq
  !
  if(mode.eq.1) then
    !
    nq=nqseg*nqbnd+1
    allocate(qvec(1:nq, 1:3))
    !
    do iq=1, nq
      get_qvec(qvec(iq), iq)
    enddo
    !
  else
    !
    
  endif
