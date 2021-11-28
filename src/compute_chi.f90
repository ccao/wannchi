SUBROUTINE compute_chi_bare_diag(chi, qv)
  !
  use constants, only : dp, eps9, cmplx_0, cmplx_i, stdout
  use para,      only : first_k, last_k, para_merge0
  use banddata,  only : kvec, nbnd, occ, eig, egv, nkpt
  use input,     only : omega, eps
  !
  implicit none
  !
  real(dp) qv(1:3)
  complex(dp) chi, chi_tmp
  !
  real(dp) ikv(1:3), jkv(1:3)
  !
  integer io, jo, ik, jk
  integer ii, jj
  !
  integer kpt_index
  !
  do ik=first_k, last_k
    ikv(:)=kvec(:, ik)
    jkv(:)=ikv(:)+qv(:)
    jk=kpt_index(jkv)
    !
    do io=1, nbnd
      do jo=1, nbnd
        !
        if (abs(occ(io, ik)-occ(jo, jk))>eps9) then
          !
          chi_tmp = cmplx_0
          do ii=1, nbnd
            do jj=1, nbnd
              chi_tmp = chi_tmp+egv(ii, io, ik)*CONJG(egv(jj, io, ik))*egv(jj, jo, jk)*CONJG(egv(ii, jo, jk))
            enddo
          enddo ! ii
          !
          chi=chi-chi_tmp*(occ(jo, jk)-occ(io, ik))/(omega+eig(jo, jk)-eig(io, ik)+cmplx_i*eps)
          !
        endif
        !
      enddo ! jo
    enddo ! io
    !
  enddo ! ik
  !
  chi=chi/nkpt
  !
  CALL para_merge0(chi)
  !
END SUBROUTINE

