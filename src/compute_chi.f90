SUBROUTINE compute_chi_diag(chi, qv)
  !
  use constants, only: dp, cmplx_0, cmplx_i
  use banddata,  only: eig, egv, nkx, nky, nkz, nkpt, nbnd, ef
  !
  implicit none
  !
  complex(dp) chi
  !
  real(dp) qv(1:3)
  !
  !
  real(dp) ikv(1:3), jkv(1:3)
  integer, external :: norm
  !
  integer ik, jk
  integer io, jo
  complex(dp) fact
  real(dp) eta
  !
  eta=1.0d-3
  !
  chi=cmplx_0
  !
  do ik=1, nkpt
    CALL map_to_kpt(ikv, ik, nkx, nky, nkz)
    jkv(:)=ikv(:)+qv(:)
    CALL map_to_idx(jk, jkv, nkx, nky, nkz)
    !
    do io=1, nbnd
      do jo=1, nbnd
        if ( (eig(io, ik)>ef).and.(eig(jo, jk)<ef) ) then
          fact=SUM(egv(:, io, ik)*conjg(egv(:, jo, jk)))
          chi=chi+conjg(fact)*fact/(eig(jo, jk)-eig(io, ik)+eta*cmplx_i)
        endif
        if ( (eig(io, ik)<ef).and.(eig(jo, jk)>ef) ) then
          fact=SUM(egv(:, io, ik)*conjg(egv(:, jo, jk)))
          chi=chi-conjg(fact)*fact/(eig(jo, jk)-eig(io, ik)-eta*cmplx_i)
        endif
      enddo
    enddo
    !
  enddo
  !
  chi=-chi/nkpt
  !
END SUBROUTINE

