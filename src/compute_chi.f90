MODULE chi_internal
  !
  use constants, only : dp, twopi, cmplx_i, cmplx_1, cmplx_0
  !
  implicit none
  !
  ! THESE VARIABLES ARE USED ONLY IN RPA CASE
  integer ndimf
  ! matrix dimension for FF-block
  integer ndimc
  ! matrix dimension for CC-block
  integer, dimension(:, :), allocatable :: idx_list_ff
  ! indice pair mapping between F and ff in full matrix
  integer, dimension(:), allocatable :: idx_list_cc
  ! indice mapping in CC block
  !complex(dp), dimension(:, :), allocatable :: USff, UCff
  !! U matrix in FF block
  !! US(:, :) is spin; UD(:, :) is charge
  !
  integer nUcp
  ! Number of U matrix elements in compact form
  real(dp), dimension(:), allocatable :: USff_cp, UCff_cp
  ! U matrix in compact form
  integer, dimension(:, :), allocatable :: idxUcp
  ! nonzero element indices in full matrix for U matrix
  ! U_{FF}(idxUcp(1,ii), idxUcp(2,ii))=UFF_{cp}(ii)
  !
  ! THESE VARIABLES ARE USED IN CHI0 CALCULATIONS
  !
  real(dp), dimension(:), allocatable :: ek, ekq, occ_k, occ_kq
  ! Eigen value and occupation at k and k+q
  real(dp), dimension(:, :), allocatable :: Skq
  ! Form factor between k & k+q
  !   Skq(ibnd, jbnd)
  !   ibnd is band index for k; jbnd is band index for kq
  complex(dp), dimension(:, :), allocatable :: hk, hkq
  ! Hamiltonian at k and k+q
  !
  contains
  !
  subroutine finalize_chi_internal()
    !
    implicit none
    !
    if (allocated(UCff_cp)) deallocate(Ucff_cp)
    if (allocated(USff_cp)) deallocate(USff_cp)
    if (allocated(idxUcp)) deallocate(idxUcp)
    if (allocated(idx_list_ff)) deallocate(idx_list_ff)
    if (allocated(idx_list_cc)) deallocate(idx_list_cc)
    !if (allocated(USff)) deallocate(USff)
    !if (allocated(UCff)) deallocate(UCff)
    if (allocated(ek)) deallocate(ek)
    if (allocated(ekq)) deallocate(ekq)
    if (allocated(occ_k)) deallocate(occ_k)
    if (allocated(occ_kq)) deallocate(occ_kq)
    if (allocated(Skq)) deallocate(Skq)
    if (allocated(hk)) deallocate(hk)
    if (allocated(hkq)) deallocate(hkq)
    nUcp=0
    ndimf=0
    ndimc=0
    !
  end subroutine

  subroutine init_chi_internal()
    !
    use constants, only : dp, stdout
    use wanndata,  only : norb
    use impurity,  only : fulldim, nimp, ndim, basis_map, nnlow, U_imp, J_imp, V_imp, Jp_imp
    use para,      only : inode
    use input,     only : level
    !
    implicit none
    !
    integer ii, jj, kk, nn, i1, i2, j1, j2, nn1, nn2, iucp
    integer, dimension(norb) :: basis_partition
    !
    allocate(ek(norb), ekq(norb), occ_k(norb), occ_kq(norb))
    allocate(Skq(norb, norb), hk(norb, norb), hkq(norb, norb))
    !
    if (level==1) then
      ndimc=norb-fulldim
      ndimf=0
      nUcp=0
      do ii=1, nimp
        ndimf=ndimf+ndim(ii)*ndim(ii)
        nUcp=nUcp+3*ndim(ii)*ndim(ii)-2*ndim(ii)
      enddo
      !
      allocate(idx_list_ff(2, ndimf), idx_list_cc(ndimc))
      !allocate(USff(ndimf, ndimf), UCff(ndimf, ndimf))
      allocate(USff_cp(nUcp), UCff_cp(nUcp), idxUcp(2, nUcp))
      !
      nn=1
      iucp=1
      do kk=1, nimp
        nn1=nn
        do ii=1, ndim(kk)
          do jj=1, ndim(kk)
            idx_list_ff(1, nn)=basis_map(nnlow(kk)+ii-1)
            idx_list_ff(2, nn)=basis_map(nnlow(kk)+jj-1)
            nn=nn+1
          enddo
        enddo
        nn2=nn-1
        !
        do ii=nn1, nn2
          i1=idx_list_ff(1, ii)
          i2=idx_list_ff(2, ii)
          do jj=nn1, nn2
            j1=idx_list_ff(1, jj)
            j2=idx_list_ff(2, jj)
            if ((i1.eq.i2).and.(j1.eq.j2)) then
              if (i1.eq.j1) then
                !USff(ii, jj)=U_imp(kk)
                !UCff(ii, jj)=U_imp(kk)
                USff_cp(iucp)=U_imp(kk)
                UCff_cp(iucp)=U_imp(kk)
                idxUcp(1, iucp)=ii
                idxUcp(2, iucp)=jj
                iucp=iucp+1
              else
                !USff(ii, jj)=0.5d0*J_imp(kk)
                !UCff(ii, jj)=2.d0*V_imp(kk)
                USff_cp(iucp)=0.5d0*J_imp(kk)
                UCff_cp(iucp)=2.d0*V_imp(kk)
                idxUcp(1, iucp)=ii
                idxUcp(2, iucp)=jj
                iucp=iucp+1
              endif
            endif
            if ((i1.eq.j1).and.(i2.eq.j2).and.(i1.ne.i2)) then
              !USff(ii, jj)=0.25d0*J_imp(kk)+V_imp(kk)
              !UCff(ii, jj)=0.75d0*J_imp(kk)-V_imp(kk)
              USff_cp(iucp)=0.25d0*J_imp(kk)+V_imp(kk)
              UCff_cp(iucp)=0.75d0*J_imp(kk)-V_imp(kk)
              idxUcp(1, iucp)=ii
              idxUcp(2, iucp)=jj
              iucp=iucp+1
            endif
            if ((i1.eq.j2).and.(i2.eq.j1).and.(i1.ne.i2)) then
              !USff(ii, jj)=Jp_imp(kk)
              !UCff(ii, jj)=Jp_imp(kk)
              USff_cp(iucp)=Jp_imp(kk)
              UCff_cp(iucp)=Jp_imp(kk)
              idxUcp(1, iucp)=ii
              idxUcp(2, iucp)=jj
              iucp=iucp+1
            endif
          enddo
        enddo
        !
      enddo ! kk
      !
      basis_partition(:)=0
      basis_partition(basis_map(:))=1
      !
      nn=1
      do ii=1, norb
        if (basis_partition(ii).eq.0) then
          idx_list_cc(nn)=ii
          nn=nn+1
        endif
      enddo
      !
      if (inode.eq.0) then
        write(stdout, *) "  FF-block indice-pair:"
        do ii=1, ndimf
          write(stdout, '(1I5,A,1I4,A,1I4,A)') ii, ' <==> (', idx_list_ff(1,ii), ',', idx_list_ff(2, ii), ')'
        enddo
        write(stdout, *) "  C-block indices:"
        do ii=1, ndimc
          write(stdout, '(1I5,A,1I5)') ii, ' <==> ', idx_list_cc(ii)
        enddo
        !
        write(stdout, *) "  Nonzero U elements:"
        write(stdout, *) " I1   I2   I3   I4       U^S        U^C"
        do ii=1, nUcp
          write(*, '(4I5,2F9.4)')  & 
      idx_list_ff(1, idxUcp(1, ii)), idx_list_ff(2, idxUcp(1, ii)), & 
      idx_list_ff(1, idxUcp(2, ii)), idx_list_ff(2, idxUcp(2, ii)), &
      USff_cp(ii), UCff_cp(ii)
        enddo
      endif
    endif ! level=1  
    !
  end subroutine
  !
END MODULE

SUBROUTINE calc_chi_bare_trace(chi0, nu, nnu, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   of frequency nu
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, beta
  use para,         only : inode, para_merge_cmplx, distribute_calc, first_idx, last_idx
  use chi_internal, only : hk, hkq
  !
  implicit none
  !
  integer nnu
  complex(dp), dimension(nnu) :: chi0, nu
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  !complex(dp), dimension(norb, norb) :: hk, hkq
  !
  integer nkpt
  integer ikx, iky, ikz, ikpt, ii
  integer inu
  !
  chi0(:)=cmplx_0
  !
  if (SUM(qv(:)*qv(:))<eps6) then
    qv(:)=0.001d0
  endif
  !
  nkpt=nkx*nky*nkz
  !
  call distribute_calc(nkpt)
  !
  do ikpt=first_idx, last_idx
    if (inode.eq.0) then
      write(stdout, '(A,1I5,A,1I5)') '   Working on Kpt: ', ikpt, ' out of ', last_idx-first_idx+1
    endif
    ikz=mod(ikpt-1, nkz)
    iky=mod((ikpt-ikz-1)/nkz, nky)
    ikx=(ikpt-1)/(nky*nkz)
    k(1)=ikx*1.d0/nkx
    k(2)=iky*1.d0/nky
    k(3)=ikz*1.d0/nkz
    kq=k+qv
    !
    call calc_hk(hk, k)
    call calc_hk(hkq, kq)
    !
    call calc_eig_occ_Skq(beta)
    !
    do inu=1, nnu
      call calc_partial_chi_bare_trace(chi0(inu), nu(inu))
    enddo
    !
  enddo
  !
  call para_merge_cmplx(chi0, nnu)
  !
  chi0(:)=chi0(:)/nkpt
  !
END SUBROUTINE

SUBROUTINE calc_partial_chi_bare_trace(chi, w)
  !
  ! This subroutine calculates part of the susceptibility,
  !   a summation over k (integration over BZ) is required to obtain full chi
  !
  use constants, only : dp, eps9, cmplx_0
  use wanndata,  only : norb
  use chi_internal, only: ek, ekq, occ_k, occ_kq, Skq
  !
  implicit none
  !
  complex(dp) :: chi, w
  ! chi: output susceptibility
  !   w: input frequency
  !
  integer ibnd, jbnd, ii, jj
  real(dp) :: occ_diff
 !
  do ibnd=1, norb
    do jbnd=1, norb
      !
      occ_diff=occ_k(ibnd)-occ_kq(jbnd)
      if (abs(occ_diff)>eps9) then
        !
        chi=chi+Skq(ibnd, jbnd)*occ_diff/(w+ekq(jbnd)-ek(ibnd))
        !
      endif
      !
    enddo ! jbnd
  enddo ! ibnd
  !
END SUBROUTINE

SUBROUTINE calc_chi_rpa_from_bare(chiRPAff, chiRPAcc, chiRPAfc, chiRPAcf, chi0ff, chi0cc, chi0fc, chi0cf, mode)
  !
  use constants, only :  dp, cmplx_0, cmplx_1, eps6
  use chi_internal, only : ndimf, ndimc, USff_cp, UCff_cp, idxUcp, nUcp
  use linalgwrap, only : invmat, matmulsparse, sparsemulmat
  !
  implicit none
  !
  integer mode
  ! mode=0: charge
  ! mode=1: spin
  complex(dp), dimension(ndimf, ndimf) :: chiRPAff, chi0ff
  complex(dp), dimension(ndimc, ndimc) :: chiRPAcc, chi0cc
  complex(dp), dimension(ndimf, ndimc) :: chiRPAfc, chi0fc
  complex(dp), dimension(ndimc, ndimf) :: chiRPAcf, chi0cf
  !
  complex(dp), dimension(ndimf, ndimf) :: Dff, Vff
  complex(dp), dimension(ndimc, ndimf) :: tmpcf
  !
  integer ii
  !
  Dff=cmplx_0
  do ii=1, ndimf
    Dff(ii, ii)=cmplx_1
  enddo
  !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  !call zgemm('N', 'N', ndimf, ndimf, ndimf, -cmplx_1, chi0ff, ndimf, USff, ndimf, cmplx_1, Dff, ndimf)
  if (mod(mode, 10).eq.1) then
    call matmulsparse(Dff, chi0ff, USff_cp, idxUcp, ndimf, nUcp, -1.d0, 1.d0)
    call invmat(Dff, ndimf)
    call sparsemulmat(Vff, USff_cp, Dff, idxUcp, ndimf, nUcp, 1.d0, 0.d0)
    call zgemm('N', 'N', ndimc, ndimf, ndimf, cmplx_1, chi0cf, ndimc, Vff, ndimf, cmplx_0, tmpcf, ndimc)
  else
    call matmulsparse(Dff, chi0ff, UCff_cp, idxUcp, ndimf, nUcp,  1.d0, 1.d0)
    call invmat(Dff, ndimf)
    call sparsemulmat(Vff, UCff_cp, Dff, idxUcp, ndimf, nUcp, 1.d0, 0.d0)
    call zgemm('N', 'N', ndimc, ndimf, ndimf, -cmplx_1, chi0cf, ndimc, Vff, ndimf, cmplx_0, tmpcf, ndimc)
  endif
  ! D^S_{FF}=(1-\chi^0_{FF}*U^S_{FF})^{-1}
  ! D^C_{FF}=(1+\chi^0_{FF}*U^C_{FF})^{-1}
  ! Vff=Uff*Dff
  ! For spin part: tmpcf=chi0cf*Vff
  ! For charge part: tmpcf=-chi0cf*Vff
  call zgemm('N', 'N', ndimf, ndimf, ndimf, cmplx_1, Dff, ndimf, chi0ff, ndimf, cmplx_0, chiRPAff, ndimf)
  ! chiff=Dff*chi0ff
  call zgemm('N', 'N', ndimf, ndimc, ndimf, cmplx_1, Dff, ndimf, chi0fc, ndimf, cmplx_0, chiRPAfc, ndimf)
  ! chifc=Dff*chi0fc
  chiRPAcf(:, :)=chi0cf(:, :)
  call zgemm('N', 'N', ndimc, ndimf, ndimf, cmplx_1, tmpcf, ndimc, chi0ff, ndimf, cmplx_1, chiRPAcf, ndimc)
  ! chiScf=chi0cf+tmpcf*chi0ff
  chiRPAcc(:, :)=chi0cc(:, :)
  call zgemm('N', 'N', ndimc, ndimc, ndimf, cmplx_1, tmpcf, ndimc, chi0fc, ndimf, cmplx_1, chiRPAcc, ndimc)
  ! chiScc=chi0cc+tmpcf*chi0fc
end subroutine

subroutine calc_chi_rpa_trace_from_matrix(chi, chiff, chicc, chifc, chicf)
  !
  use constants, only : dp, cmplx_0
  use chi_internal, only : ndimf, ndimc, idx_list_ff
  !
  implicit none
  !
  complex(dp) :: chi
  complex(dp), dimension(ndimf, ndimf) :: chiff
  complex(dp), dimension(ndimf, ndimc) :: chifc
  complex(dp), dimension(ndimc, ndimc) :: chicc
  complex(dp), dimension(ndimc, ndimf) :: chicf
  !
  integer ii, jj
  !
  chi=cmplx_0
  !
  do ii=1, ndimf
    if (idx_list_ff(1, ii).eq.idx_list_ff(2, ii)) then
      do jj=1, ndimf
        if (idx_list_ff(1, jj).eq.idx_list_ff(2, jj)) then
          chi=chi+chiff(ii, jj)
        endif
      enddo
      do jj=1, ndimc
        chi=chi+chifc(ii, jj)+chicf(jj, ii)
      enddo
    endif
  enddo
  do ii=1, ndimc
    do jj=1, ndimc
      chi=chi+chicc(ii, jj)
    enddo
  enddo
  !
end subroutine

SUBROUTINE calc_chi_rpa_trace(chis, chic, chi0, w, nnu, qv)
  !
  use constants,    only : dp, cmplx_0, cmplx_1, eps6
  use wanndata,     only : norb
  use chi_internal, only : ndimf, ndimc, idx_list_cc, idx_list_ff, hk, hkq
  use para,         only : distribute_calc, inode, first_idx, last_idx, para_merge_cmplx
  use input,        only : beta, nkx, nky, nkz
  !
  implicit none
  !
  integer nnu
  complex(dp), dimension(nnu), intent(out) :: chis, chic, chi0
  complex(dp), dimension(nnu), intent(in) :: w
  real(dp), dimension(3)  :: qv
  complex(dp), dimension(ndimf, ndimf, nnu) :: chi0ff
  complex(dp), dimension(ndimc, ndimc, nnu) :: chi0cc
  complex(dp), dimension(ndimf, ndimc, nnu) :: chi0fc
  complex(dp), dimension(ndimc, ndimf, nnu) :: chi0cf
  complex(dp), dimension(ndimf, ndimf) :: chiSff, chiCff
  complex(dp), dimension(ndimc, ndimc) :: chiScc, chiCcc
  complex(dp), dimension(ndimf, ndimc) :: chiSfc, chiCfc
  complex(dp), dimension(ndimc, ndimf) :: chiScf, chiCcf
  !
  complex(dp), dimension(norb, norb) :: Uk, Ukq
  !
  integer ikpt, nkpt, ikx, iky, ikz, ii, jj, inu
  real(dp), dimension(3) :: k, kq
  !
  if (SUM(qv(:)*qv(:))<eps6) then
    qv(:)=0.001d0
  endif
  !
  nkpt=nkx*nky*nkz
  !
  call distribute_calc(nkpt)
  !
  chi0ff=cmplx_0
  chi0cc=cmplx_0
  chi0cf=cmplx_0
  chi0fc=cmplx_0
  !
  do ikpt=first_idx, last_idx
    ikz=mod(ikpt-1, nkz)
    iky=mod((ikpt-ikz-1)/nkz, nky)
    ikx=(ikpt-1)/(nky*nkz)
    k(1)=ikx*1.d0/nkx
    k(2)=iky*1.d0/nky
    k(3)=ikz*1.d0/nkz
    kq=k+qv
    !
    call calc_hk(hk, k)
    call calc_hk(hkq, kq)
    !
    call calc_eig_occ_Umat(Uk, Ukq, beta)
    !
    do inu=1, nnu
      call calc_partial_chi0_matrix(chi0ff(:, :, inu), chi0cc(:, :, inu), chi0fc(:, :, inu), chi0cf(:, :, inu), w(inu), Uk, Ukq)
      !call calc_partial_chi0_matrix(chi0ff_part, chi0cc_part, chi0fc_part, chi0cf_part, w, Uk, Ukq)
    enddo
  enddo
  !
  call para_merge_cmplx(chi0ff, ndimf*ndimf*nnu)
  call para_merge_cmplx(chi0cc, ndimc*ndimc*nnu)
  call para_merge_cmplx(chi0fc, ndimf*ndimc*nnu)
  call para_merge_cmplx(chi0cf, ndimc*ndimf*nnu)
  !
  chi0ff=chi0ff/nkpt
  chi0cc=chi0cc/nkpt
  chi0fc=chi0fc/nkpt
  chi0cf=chi0cf/nkpt
  !
  chis=cmplx_0
  chic=cmplx_0
  chi0=cmplx_0
  !
  call distribute_calc(nnu)
  !
  do inu=first_idx, last_idx
    !
    call calc_chi_rpa_from_bare(chiSff, chiScc, chiSfc, chiScf, chi0ff(:, :, inu), chi0cc(:, :, inu), chi0fc(:, :, inu), chi0cf(:, :, inu), 1)
    call calc_chi_rpa_from_bare(chiCff, chiCcc, chiCfc, chiCcf, chi0ff(:, :, inu), chi0cc(:, :, inu), chi0fc(:, :, inu), chi0cf(:, :, inu), 0)
    !
    ! Trace
    !
    call calc_chi_rpa_trace_from_matrix(chi0(inu), chi0ff(:, :, inu), chi0cc(:, :, inu), chi0fc(:, :, inu), chi0cf(:, :, inu))
    call calc_chi_rpa_trace_from_matrix(chiS(inu), chiSff, chiScc, chiSfc, chiScf)
    call calc_chi_rpa_trace_from_matrix(chiC(inu), chiCff, chiCcc, chiCfc, chiCcf)
    !
  enddo
  !
  call para_merge_cmplx(chi0, nnu)
  call para_merge_cmplx(chiS, nnu)
  call para_merge_cmplx(chiC, nnu)
  !
END SUBROUTINE

SUBROUTINE calc_partial_chi0_matrix(chi0ff, chi0cc, chi0fc, chi0cf, w, Uk, Ukq)
  !
  use constants, only : dp, cmplx_0, cmplx_1, eps9
  use wanndata,  only : norb
  use chi_internal, only : ndimf, ndimc, idx_list_cc, idx_list_ff, ek, ekq, occ_k, occ_kq
  !
  implicit none
  !  
  complex(dp), dimension(ndimf, ndimf), intent(out) :: chi0ff
  complex(dp), dimension(ndimc, ndimc), intent(out) :: chi0cc
  complex(dp), dimension(ndimf, ndimc), intent(out) :: chi0fc
  complex(dp), dimension(ndimc, ndimf), intent(out) :: chi0cf
  complex(dp), dimension(norb, norb), intent(in) :: Uk, Ukq
  complex(dp), intent(in) :: w
  !
  integer ii, jj, i1, i2, j1, j2
  integer iorb, jorb
  real(dp) :: occ_diff
  complex(dp) :: fact
  !
  !chi0ff=cmplx_0
  !chi0cc=cmplx_0
  !chi0fc=cmplx_0
  !chi0cf=cmplx_0
  !
  ! Obtain the Full matrix of chi0_ff
  !
  do iorb=1, norb   ! mu
    do jorb=1, norb ! nu
      occ_diff=occ_k(iorb)-occ_kq(jorb)
      if (abs(occ_diff)>eps9) then
        fact=occ_diff/(w+ekq(jorb)-ek(iorb))
        do ii=1, ndimf ! upper indices
          i1=idx_list_ff(1, ii) !p
          i2=idx_list_ff(2, ii) !q
          do jj=1, ndimf ! lower indices
            j1=idx_list_ff(1, jj) ! s
            j2=idx_list_ff(2, jj) ! t
            chi0ff(ii, jj)=chi0ff(ii, jj)+fact*Uk(j1, iorb)*conjg(Uk(i1, iorb))*Ukq(i2, jorb)*conjg(Ukq(j2, jorb))
          enddo
        enddo
        !
        do ii=1, ndimf ! upper indices
          i1=idx_list_ff(1, ii) !p
          i2=idx_list_ff(2, ii) !q
          do jj=1, ndimc
            j1=idx_list_cc(jj)  !s,t
            chi0fc(ii, jj)=chi0fc(ii, jj)+fact*Uk(j1, iorb)*conjg(Uk(i1, iorb))*Ukq(i2, jorb)*conjg(Ukq(j1, jorb))
          enddo
        enddo
        !
        do ii=1, ndimc ! 
          i1=idx_list_cc(ii) ! p,q
          do jj=1, ndimf
            j1=idx_list_ff(1, jj) ! s
            j2=idx_list_ff(2, jj) ! t
            chi0cf(ii, jj)=chi0cf(ii, jj)+fact*Uk(j1, iorb)*conjg(Uk(i1, iorb))*Ukq(i1, jorb)*conjg(Ukq(j2, jorb))
          enddo
        enddo
        !
        do ii=1, ndimc !
          i1=idx_list_cc(ii) ! p,q
          do jj=1, ndimc !
            j1=idx_list_cc(jj) ! s,t
            chi0cc(ii, jj)=chi0cc(ii, jj)+fact*Uk(j1, iorb)*conjg(Uk(i1, iorb))*Ukq(i1, jorb)*conjg(Ukq(j1, jorb))
          enddo
        enddo
      endif
    enddo ! jorb
  enddo! iorb
  !
END SUBROUTINE

SUBROUTINE calc_chi_corr_trace(chi, chi0, nnu, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   Up to first nnu Matsubara frequencies
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, nexact, beta
  use para,         only : inode
  use chi_internal, only : hk, hkq
  !
  implicit none
  !
  integer      :: nnu
  complex(dp), dimension(nnu) :: chi, chi0
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  !
  complex(dp), dimension(nnu) :: chi_part, chi0_part
  !
  integer ikx, iky, ikz
  !
  chi(:)=cmplx_0
  chi0(:)=cmplx_0
  !
  if (SUM(qv(:)*qv(:))<eps6) then
    qv(:)=0.001d0
  endif
  !
  do ikx=0, nkx-1
    k(1)=ikx*1.d0/nkx
    do iky=0, nky-1
      k(2)=iky*1.d0/nky
      do ikz=0, nkz-1
        k(3)=ikz*1.d0/nkz
        kq=k+qv
        !
        call calc_hk(hk, k)
        call calc_hk(hkq, kq)
        !
        call calc_eig_occ_Skq(beta)
        !
        call calc_partial_chi_corr_trace(chi_part, chi0_part, nnu, nexact, beta)
        !
        chi(:)=chi(:)+chi_part(:)
        chi0(:)=chi0(:)+chi0_part(:)
        !
      enddo ! ikz
    enddo
  enddo
  !
  chi(:)=chi(:)/(nkx*nky*nkz)
  chi0(:)=chi0(:)/(nkx*nky*nkz)
  !
END SUBROUTINE

SUBROUTINE calc_partial_chi_corr_trace(chi, chi0, nnu, nexact, beta)
  !
  ! This subroutine calculates part of the susceptibility (one kpt in the K-mesh),
  !   Parallel over frequency
  !
  use constants, only : dp, eps9, cmplx_0, twopi, cmplx_i, stdout
  use wanndata,  only : norb
  use para,      only : distribute_calc, first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx
  use chi_internal, only : ek, ekq, Skq, hk, hkq
  !
  implicit none
  !
  integer, intent(in)       :: nnu
  ! First nnu Matsubara frequencies will be calculated
  integer, intent(in)       :: nexact
  ! We only calculate -nexact to nexact contributions exactly (the difference decays very fast)
  real(dp), intent(in)      :: beta
  complex(dp), dimension(nnu), intent(out) :: chi, chi0
  ! chi: correlated susceptibility
  ! chi0: noninteracting susceptibility
  !
  complex(dp) :: w1, w2, nu
  integer ibnd, jbnd
  integer iom, inu
  ! In Matsubara calculations, displacement of frequency points
  complex(dp), dimension(:, :, :), allocatable :: gk, gkq
  ! Interacting Green's function (local, for calculations)
  complex(dp), dimension(:, :, :), allocatable :: gk_full, gkq_full
  ! Interacting Green's function (complete, only relevant on master node)
  complex(dp)                 :: chi_diff
  !
  ! Calculate chi0 first
  !
  do inu=0, nnu-1
    nu=inu*twopi/beta*cmplx_i
    call calc_partial_chi_bare_trace(chi0(inu+1), nu)
  enddo
  !
  ! Prepare calculation of correlated chi
  !   gk_full, gkq_full contains all frequencies
  if (inode.eq.0) then
    !allocate(gk_full(norb, norb, 2*nfreq), gkq_full(norb, norb, 2*nfreq))
    allocate(gk_full(norb, norb, 2*nexact), gkq_full(norb, norb, 2*nexact))
  else
    allocate(gk_full(1,1,1), gkq_full(1,1,1))
  endif
  !
  !call distribute_calc(2*nfreq)
  call distribute_calc(2*nexact)
  !
  allocate(gk(norb, norb, last_idx-first_idx+1))
  allocate(gkq(norb, norb, last_idx-first_idx+1))
  !
  do iom=first_idx, last_idx
    !call calc_corr_matsgf(gk(:, :, iom-first_idx+1), hk, iom-nfreq, .false.)
    !call calc_corr_matsgf(gkq(:, :, iom-first_idx+1), hkq, iom-nfreq, .false.)
    call calc_corr_matsgf(gk(:, :, iom-first_idx+1), hk, iom-nexact, .false.)
    call calc_corr_matsgf(gkq(:, :, iom-first_idx+1), hkq, iom-nexact, .false.)
  enddo
  !
  call para_collect_cmplx(gk_full, gk, norb*norb)
  call para_collect_cmplx(gkq_full, gkq, norb*norb)
  !
  deallocate(gk, gkq) ! Free the memory
  !
  ! Up till now, gk_full and gkq_full should contain all green's functions
  ! Now do the actual calculation
  !
  chi(:)=cmplx_0
  do inu=0, nnu-1
    !
    !nu=inu*twopi/beta*cmplx_i
    if (inode.eq.0) write(stdout, '(A,1I8)') "    Done inu :", inu
    !
    !call distribute_calc(2*nfreq-inu)
    call distribute_calc(2*nexact-inu)
    !
    allocate(gk(norb, norb, last_idx-first_idx+1))
    allocate(gkq(norb, norb, last_idx-first_idx+1))
    !
    !call para_distribute_cmplx(gk_full(:, :, 1:2*nfreq-inu), gk, norb*norb)
    !call para_distribute_cmplx(gkq_full(:, :, 1+inu:2*nfreq), gkq, norb*norb)
    call para_distribute_cmplx(gk_full(:, :, 1:2*nexact-inu), gk, norb*norb)
    call para_distribute_cmplx(gkq_full(:, :, 1+inu:2*nexact), gkq, norb*norb)
    !
    do iom=first_idx, last_idx
      !
      !w1=(iom-nfreq-0.5d0)*twopi/beta*cmplx_i
      !w2=(iom-nfreq+inu-0.5d0)*twopi/beta*cmplx_i
      w1=(iom-nexact-0.5d0)*twopi/beta*cmplx_i
      w2=(iom-nexact+inu-0.5d0)*twopi/beta*cmplx_i
      !
      do ibnd=1, norb
        do jbnd=1, norb
          chi_diff=gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)-Skq(ibnd, jbnd)/((w1-ek(ibnd))*(w2-ekq(jbnd)))
          chi(inu+1)=chi(inu+1)-chi_diff
        enddo
      enddo ! ibnd
      !
    enddo ! iom
    !
    deallocate(gk, gkq)
    !
  enddo ! inu
  !
  deallocate(gk_full, gkq_full)
  !
  call para_merge_cmplx(chi, nnu)
  chi(:)=chi(:)/beta+chi0(:)
  if (inode.eq.0) write(stdout, *)
  !
END SUBROUTINE

SUBROUTINE calc_eig_occ_Skq(beta)
  !
  !
  use constants, only : dp, cmplx_1, cmplx_0
  use wanndata,  only : norb
  use linalgwrap, only : eigen
  use chi_internal, only : ek, ekq, occ_k, occ_kq, Skq, hk, hkq
  !
  implicit none
  !
  real(dp), intent(in) :: beta
  !
  complex(dp), dimension(norb, norb) :: Uk, Ukq, Stmp
  real(dp) calc_occ
  !
  integer ii, jj
  !
  Uk=hk
  call eigen(ek, Uk, norb)
  !
  Ukq=hkq
  call eigen(ekq, Ukq, norb)
  !
  do ii=1, norb
    occ_k(ii)=calc_occ(ek(ii), beta)
    occ_kq(ii)=calc_occ(ekq(ii), beta)
  enddo
  !
  ! calculates the structural function, 
  !  Skq(ibnd, jbnd)= | \sum_s conjg(egv_k(s, ibnd))*egv_kq(s, jbnd) |^2
  !   or in matrix form: Skq = | U_k^dagger * U_kq |^2
  !   !!! PLEASE BE AWARE THAT ibnd/jbnd order is relevant
  !       ibnd is the band index for k
  !    and jbnd is the band index for kq
  !
  !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call zgemm('C', 'N', norb, norb, norb, cmplx_1, Uk, norb, Ukq, norb, cmplx_0, Stmp, norb)
  !
  do ii=1, norb
    do jj=1, norb
      !Stmp(ii, jj)=sum(conjg(Uk(:, ii))*Ukq(:, jj))
      Skq(ii, jj)=abs(Stmp(ii, jj))**2
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE analyze_part_chi_corr_trace(chi, chi0, nnu, beta)
  !
  !  This subroutine analyzes the composition of chi_corr
  !    Please be noted that unlike the calculation subroutine
  !      nnu here is the matsubara frequency being analyzed
  !      and chi/chi0 contains the frequency components instead of final results
  !
  use constants, only : dp, cmplx_0, twopi, cmplx_i
  use wanndata,  only : norb
  use impurity,  only : nfreq
  use para,      only : first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx, distribute_calc
  use chi_internal, only : ek, ekq, occ_k, occ_kq, hk, hkq, Skq
  !
  implicit none
  !
  integer       :: nnu
  complex(dp), dimension(2*nfreq-nnu) :: chi, chi0
  real(dp)      :: beta
  !
  complex(dp)  :: w1, w2
  integer  ibnd, jbnd
  integer  iom
  !
  complex(dp), dimension(:, :, :), allocatable :: gk, gkq
  complex(dp), dimension(:, :, :), allocatable :: gk_full, gkq_full
  !
  chi=cmplx_0
  chi0=cmplx_0
  !
  if (inode.eq.0) then
    allocate(gk_full(norb, norb, 2*nfreq), gkq_full(norb, norb, 2*nfreq))
  else
    allocate(gk_full(1,1,1), gkq_full(1,1,1))
  endif
  !
  call distribute_calc(2*nfreq)
  allocate(gk(norb, norb, last_idx-first_idx+1))
  allocate(gkq(norb, norb, last_idx-first_idx+1))
  !
  do iom=first_idx, last_idx
    call calc_corr_matsgf(gk(:, :, iom-first_idx+1), hk, iom-nfreq, .false.)
    call calc_corr_matsgf(gkq(:, :, iom-first_idx+1), hkq, iom-nfreq, .false.)
  enddo
  !
  call para_collect_cmplx(gk_full, gk, norb*norb)
  call para_collect_cmplx(gkq_full, gkq, norb*norb)
  !
  deallocate(gk, gkq) ! Free the memory
  !
  call distribute_calc(2*nfreq-nnu)
  !
  allocate(gk(norb, norb, last_idx-first_idx+1))
  allocate(gkq(norb, norb, last_idx-first_idx+1))
  !
  call para_distribute_cmplx(gk_full(:, :, 1:2*nfreq-nnu), gk, norb*norb)
  call para_distribute_cmplx(gkq_full(:, :, 1+nnu:2*nfreq), gkq, norb*norb)
  !
  do iom=first_idx, last_idx
    !
    w1=(iom-nfreq-0.5d0)*twopi/beta*cmplx_i
    w2=(iom-nfreq+nnu-0.5d0)*twopi/beta*cmplx_i
    !
    do ibnd=1, norb
      do jbnd=1, norb
        chi(iom)=chi(iom)+gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)
        chi0(iom)=chi0(iom)+Skq(ibnd, jbnd)/((w1-ek(ibnd))*(w2-ekq(jbnd)))
      enddo
    enddo ! ibnd
    !
  enddo ! iom
  !
  deallocate(gk, gkq)
  deallocate(gk_full, gkq_full)
  !
  call para_merge_cmplx(chi0, 2*nfreq-nnu)
  call para_merge_cmplx(chi, 2*nfreq-nnu)
  chi0=chi0/beta
  chi=chi/beta
  !
END SUBROUTINE

SUBROUTINE analyze_part_chi_corr_trace1(chi, chif, chic, nnu, beta)
  !
  !  This subroutine analyzes the composition of chi_corr
  !    Please be noted that unlike the calculation subroutine
  !      nnu here is the matsubara frequency being analyzed
  !      and chi/chi0 contains the frequency components instead of final results
  !
  use constants, only : dp, cmplx_0, twopi, cmplx_i
  use wanndata,  only : norb
  use impurity,  only : nfreq, fulldim, basis_map
  use para,      only : first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx, distribute_calc
  use chi_internal, only : hk, hkq
  !
  implicit none
  !
  integer       :: nnu
  real(dp)      :: beta
  complex(dp), dimension(2*nfreq-nnu) :: chi, chif, chic
  !
  complex(dp)  :: w1, w2
  integer  ibnd, jbnd
  integer  iom
  !
  complex(dp), dimension(:, :, :), allocatable :: gk, gkq
  complex(dp), dimension(:, :, :), allocatable :: gk_full, gkq_full
  !
  integer, dimension(norb) :: orb_type
  !
  chi=cmplx_0
  chif=cmplx_0
  chic=cmplx_0
  !
  orb_type=0
  do ibnd=1, fulldim
    orb_type(basis_map(ibnd))=1
  enddo
  !
  if (inode.eq.0) then
    allocate(gk_full(norb, norb, 2*nfreq), gkq_full(norb, norb, 2*nfreq))
  else
    allocate(gk_full(1,1,1), gkq_full(1,1,1))
  endif
  !
  call distribute_calc(2*nfreq)
  allocate(gk(norb, norb, last_idx-first_idx+1))
  allocate(gkq(norb, norb, last_idx-first_idx+1))
  !
  do iom=first_idx, last_idx
    call calc_corr_matsgf(gk(:, :, iom-first_idx+1), hk, iom-nfreq, .false.)
    call calc_corr_matsgf(gkq(:, :, iom-first_idx+1), hkq, iom-nfreq, .false.)
  enddo
  !
  call para_collect_cmplx(gk_full, gk, norb*norb)
  call para_collect_cmplx(gkq_full, gkq, norb*norb)
  !
  deallocate(gk, gkq) ! Free the memory
  !
  call distribute_calc(2*nfreq-nnu)
  !
  allocate(gk(norb, norb, last_idx-first_idx+1))
  allocate(gkq(norb, norb, last_idx-first_idx+1))
  !
  call para_distribute_cmplx(gk_full(:, :, 1:2*nfreq-nnu), gk, norb*norb)
  call para_distribute_cmplx(gkq_full(:, :, 1+nnu:2*nfreq), gkq, norb*norb)
  !
  do iom=first_idx, last_idx
    !
    w1=(iom-nfreq-0.5d0)*twopi/beta*cmplx_i
    w2=(iom-nfreq+nnu-0.5d0)*twopi/beta*cmplx_i
    !
    do ibnd=1, norb
      do jbnd=1, norb
        chi(iom)=chi(iom)+gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)
        if ((orb_type(ibnd).eq.1).and.(orb_type(jbnd).eq.1)) then
          chif(iom)=chif(iom)+gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)
        endif
        if ((orb_type(ibnd).eq.0).and.(orb_type(jbnd).eq.0)) then
          chic(iom)=chic(iom)+gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)
        endif
      enddo
    enddo ! ibnd
    !
  enddo ! iom
  !
  deallocate(gk, gkq)
  deallocate(gk_full, gkq_full)
  !
  call para_merge_cmplx(chic, 2*nfreq-nnu)
  call para_merge_cmplx(chif, 2*nfreq-nnu)
  call para_merge_cmplx(chi, 2*nfreq-nnu)
  chic=chic/beta
  chif=chif/beta
  chi=chi/beta
  !
END SUBROUTINE

SUBROUTINE analyze_chi_corr_trace(chi, chi0, nnu, qv)
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, beta
  use impurity,     only : nfreq
  use para,         only : inode
  use chi_internal, only : hk, hkq
  !
  implicit none
  !
  integer     :: nnu
  !
  complex(dp), dimension(2*nfreq-nnu) :: chi, chi0
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  complex(dp), dimension(2*nfreq) :: chi_part, chi0_part
  !
  integer ikx, iky, ikz
  !
  chi(:)=cmplx_0
  chi0(:)=cmplx_0
  !
  if (SUM(qv(:)*qv(:))<eps6) then
    qv(:)=0.001d0
  endif
  !
  do ikx=0, nkx-1
    k(1)=ikx*1.d0/nkx
    do iky=0, nky-1
      k(2)=iky*1.d0/nky
      do ikz=0, nkz-1
        k(3)=ikz*1.d0/nkz
        kq=k+qv
        !
        if (inode.eq.0) then
          write(stdout, '(A,3F9.4)') ' Working on k: ', k
        endif
        !
        call calc_hk(hk, k)
        call calc_hk(hkq, kq)
        !
        call calc_eig_occ_Skq(beta)
        !
        call analyze_part_chi_corr_trace(chi_part, chi0_part, nnu)
        !
        chi(:)=chi(:)+chi_part(:)
        chi0(:)=chi0(:)+chi0_part(:)
        !
      enddo ! ikz
    enddo
  enddo
  !
  chi(:)=chi(:)/(nkx*nky*nkz)
  chi0(:)=chi0(:)/(nkx*nky*nkz)
  !
END SUBROUTINE

SUBROUTINE analyze_chi_corr_trace1(chi, chif, chic, nnu, qv)
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, beta
  use impurity,     only : nfreq
  use para,         only : inode
  use chi_internal, only : hk, hkq
  !
  implicit none
  !
  integer, intent(in)     :: nnu
  !
  complex(dp), dimension(2*nfreq-nnu) :: chi, chif, chic
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  !
  complex(dp), dimension(2*nfreq) :: chi_part, chif_part, chic_part
  !
  integer ikx, iky, ikz
  !
  chi(:)=cmplx_0
  chif(:)=cmplx_0
  chic(:)=cmplx_0
  !
  if (SUM(qv(:)*qv(:))<eps6) then
    qv(:)=0.001d0
  endif
  !
  do ikx=0, nkx-1
    k(1)=ikx*1.d0/nkx
    do iky=0, nky-1
      k(2)=iky*1.d0/nky
      do ikz=0, nkz-1
        k(3)=ikz*1.d0/nkz
        kq=k+qv
        !
        if (inode.eq.0) then
          write(stdout, '(A,3F9.4)') ' Working on k: ', k
        endif
        !
        call calc_hk(hk, k)
        call calc_hk(hkq, kq)
        !
        call calc_eig_occ_Skq(beta)
        !
        call analyze_part_chi_corr_trace1(chi_part, chif_part, chic_part, nnu)
        !
        chi(:)=chi(:)+chi_part(:)
        chif(:)=chif(:)+chif_part(:)
        chic(:)=chic(:)+chic_part(:)
        !
      enddo ! ikz
    enddo
  enddo
  !
  chi(:)=chi(:)/(nkx*nky*nkz)
  chif(:)=chif(:)/(nkx*nky*nkz)
  chic(:)=chic(:)/(nkx*nky*nkz)
  !
END SUBROUTINE

REAL(DP) FUNCTION calc_occ(ek, beta)
  !
  use constants, only : dp
  !
  implicit none
  !
  real(dp) :: ek, beta
  !
  real(dp) :: x
  x=ek*beta
  if (x>30.d0) then
    calc_occ=0.d0
  else if (x<-30.d0) then
    calc_occ=1.d0
  else
    calc_occ=1.d0/(exp(x)+1.d0)
  endif
  !
END FUNCTION

SUBROUTINE calc_eig_occ_Umat(Uk, Ukq, beta)
  !
  !
  use constants, only : dp, cmplx_1, cmplx_0
  use wanndata,  only : norb
  use linalgwrap, only : eigen
  use chi_internal, only : ek, ekq, occ_k, occ_kq, hk, hkq
  !
  implicit none
  !
  complex(dp), dimension(norb, norb), intent(out) :: Uk, Ukq
  real(dp), intent(in) :: beta
  !
  real(dp) calc_occ
  !
  integer ii, jj
  !
  Uk=hk
  call eigen(ek, Uk, norb)
  !
  Ukq=hkq
  call eigen(ekq, Ukq, norb)
  !
  do ii=1, norb
    occ_k(ii)=calc_occ(ek(ii), beta)
    occ_kq(ii)=calc_occ(ekq(ii), beta)
  enddo
  !
END SUBROUTINE

