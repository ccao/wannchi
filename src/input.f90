MODULE input
  !
  use constants
  !
  implicit none
  !
  integer code   ! 0: wannchi 1: wanndos
  integer mode   ! Mode for calculation, 0 (default): plane; 1: band
  real(dp), allocatable :: bnd_q(:, :) ! special Q-points
  real(dp) temp  ! Temperature
  real(dp) omega ! omega ( real axis energy )
  real(dp) hubbard_u, hubbard_j, hubbard_v, hubbard_jp
  integer nqseg   ! Num of q-points per seg
  integer nqbnd   ! Num of special Q points
  !
  logical lrpa
  !
CONTAINS

SUBROUTINE finalize_input()
  !
  implicit none
  !
  if(allocated(bnd_q)) deallocate(bnd_q)
  !
END SUBROUTINE

SUBROUTINE get_qvec(qv, iq)
  !
  implicit none
  !
  integer iq
  real(dp) qv(1:3)
  integer iseg, iqb
  !
  iqb=mod(iq-1,nqbnd)
  iseg=(iq-1-iqb)/nqbnd
  !
  qv(:)=bnd_q(iseg+1,:)+iqb*(bnd_q(iseg+2,:)-bnd_q(iseg+1,:))/nqbnd
  !
END SUBROUTINE

END MODULE

