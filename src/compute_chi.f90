SUBROUTINE compute_chi_bare(qv)
  !
  use constants, only: dp, cmplx_0, cmplx_i, eps, eps4
  use banddata,  only: eig, egv, nkx, nky, nkz, nkpt, nbnd, ef
  use input,     only: temp, omega
  use chidata,   only: chi_loc
  !
  implicit none
  !
  real(dp) qv(1:3)
  !
  !
  real(dp) ikv(1:3), jkv(1:3)   !  k & k+q
  !
  real(dp) fermi_dirac
  !
  integer ii, jj, kk, ll
  integer io, jo, ik, jk
  complex(dp) fact
  real(dp) occ_io, occ_jo
  !
  if (.not.allocated(chi_loc)) allocate(chi_loc(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  !
  chi_loc(:,:,:,:)=cmplx_0
  !
  do ik=1, nkpt
    CALL map_to_kpt(ikv, ik, nkx, nky, nkz)
    jkv(:)=ikv(:)+qv(:)
    CALL map_to_idx(jk, jkv, nkx, nky, nkz)
    !
    do io=1, nbnd
      !
      occ_io=fermi_dirac(eig(io, ik)-ef, temp)
      !
      do jo=1, nbnd
        !
        occ_jo=fermi_dirac(eig(jo, jk)-ef, temp)
        !
        if (abs(occ_io-occ_jo)>eps) then
          do ii=1, nbnd
            do jj=1, nbnd
              do kk=1, nbnd
                do ll=1, nbnd
                  chi_loc(ii, jj, kk, ll) = chi_loc(ii, jj, kk, ll) - (occ_jo-occ_io)*CONJG(egv(ii, io, ik))*egv(jj, jo, jk)*egv(kk, io, ik)*CONJG(egv(ll, jo, jk))/(omega+eig(jo, jk)-eig(io, ik)+eps4*cmplx_i)/nkpt
                enddo  ! ll
              enddo  ! kk
            enddo  ! jj
          enddo  ! ii
        endif
        !
      enddo ! jo
    enddo ! io
    !
  enddo ! ik
  !
END SUBROUTINE

SUBROUTINE compute_U_mat
  !
  use constants, only: dp
  use banddata,  only: nbnd
  use input,     only: hubbard_u, hubbard_j, hubbard_v, hubbard_jp
  use chidata,   only: u_mat
  !
  implicit none
  !
  integer ii, jj
  !
  if (.not.allocated(u_mat)) allocate(u_mat(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  u_mat(:,:,:,:)=0.d0
  !
  do ii=1, nbnd
    u_mat(ii, ii, ii, ii)=hubbard_u
    do jj=1, nbnd
      u_mat(ii, ii, jj, jj)=hubbard_j/2.d0
      u_mat(ii, jj, ii, jj)=hubbard_j/4.d0+hubbard_v
      u_mat(jj, ii, ii, jj)=hubbard_jp
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE compute_chi_rpa(qv)
  !
  use para, only: inode
  use constants, only: dp, cmplx_0, cmplx_i
  use banddata,  only: nbnd
  use chidata,   only: chi_loc, chi_rpa, u_mat, u_chi, chi_tmp
  !
  implicit none
  !
  real(dp) qv(1:3)
  !
  complex(dp), allocatable :: chi_new(:, :, :, :)
  !
  integer i_loop, n_loop
  integer ii
  logical converged
  !
  n_loop = 200
  !
  if (.not.allocated(chi_rpa)) allocate(chi_rpa(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  if (.not.allocated(u_chi))  allocate(u_chi(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  allocate(chi_new(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  if (.not.allocated(chi_tmp)) allocate(chi_tmp(1:nbnd, 1:nbnd, 1:nbnd, 1:nbnd))
  !
  CALL compute_chi_bare(qv)
  !
  chi_rpa=chi_loc
  !
  CALL matrix_multiply(u_chi, u_mat, chi_loc, nbnd)
  !
  do i_loop=1, n_loop
    CALL matrix_multiply(chi_tmp, chi_rpa, u_chi, nbnd)
    chi_new=chi_loc+chi_tmp
    chi_tmp=chi_new-chi_rpa
    if (converged) exit
    chi_rpa=chi_new
  enddo
  !
  deallocate(chi_new)
  !
END SUBROUTINE

SUBROUTINE compute_chi_diag(chi, qv)
  !
  use constants, only: dp, cmplx_0
  use banddata,  only: nbnd
  use chidata,   only: chi_loc, chi_rpa
  use input,     only: lrpa
  !
  implicit none
  !
  complex(dp) chi
  !
  real(dp) qv(1:3)
  !
  integer ii, jj
  !
  chi=cmplx_0
  !
  if (lrpa) then
    CALL compute_chi_rpa(qv)
  else
    CALL compute_chi_bare(qv)
  endif
  !
  do ii=1, nbnd
    do jj=1, nbnd
      if(lrpa) then
        chi=chi+chi_rpa(ii, ii, jj, jj)
      else
        chi=chi+chi_loc(ii, ii, jj, jj)
      endif
    enddo
  enddo
  !
END SUBROUTINE

