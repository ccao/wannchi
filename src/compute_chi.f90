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
        if ( (eig(ik, io)>ef).and.(eig(jk, jo)<ef) ) then
          fact=SUM(egv(ik, :, io)*conjg(egv(jk, :, jo)))
          chi=chi+conjg(fact)*fact/(eig(jk, jo)-eig(ik, io)+eta*cmplx_i)
        else 
          if ( (eig(ik, io)<ef).and.(eig(jk, jo)>ef) ) then
            fact=SUM(egv(ik, :, io)*conjg(egv(jk, :, jo)))
            chi=chi-conjg(fact)*fact/(eig(jk, jo)-eig(ik, io)-eta*cmplx_i)
          endif
        endif
      enddo
    enddo
    !
  enddo
  !
  chi=-chi/nkpt
  !
END SUBROUTINE

