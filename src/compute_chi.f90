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

SUBROUTINE calc_eig_occ_Skq(eig_k, eig_kq, occ_k, occ_kq, Skq, hk, hkq, beta)
  !
  !
  use constants, only : dp, cmplx_1, cmplx_0
  use wanndata,  only : norb
  use linalgwrap, only : eigen
  !
  implicit none
  !
  real(dp), dimension(norb) :: eig_k, eig_kq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
  complex(dp), dimension(norb, norb) :: hk, hkq
  real(dp) :: beta
  !
  complex(dp), dimension(norb, norb) :: Uk, Ukq, Stmp
  real(dp) calc_occ
  !
  integer ii, jj
  !
  Uk=hk
  call eigen(eig_k, Uk, norb)
  !
  Ukq=hkq
  call eigen(eig_kq, Ukq, norb)
  !
  do ii=1, norb
    occ_k(ii)=calc_occ(eig_k(ii), beta)
    occ_kq(ii)=calc_occ(eig_kq(ii), beta)
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
      Skq(ii, jj)=abs(Stmp(ii, jj))**2
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_partial_chi_bare_trace(chi, w, eig_k, eig_kq, occ_k, occ_kq, Skq)
  !
  ! This subroutine calculates part of the susceptibility,
  !   a summation over k (integration over BZ) is required to obtain full chi
  !
  use constants, only : dp, eps9, cmplx_0
  use wanndata,  only : norb
  !
  implicit none
  !
  complex(dp) :: chi, w
  ! chi: output susceptibility
  !   w: input frequency
  real(dp), dimension(norb) :: eig_k, eig_kq, occ_k, occ_kq
  ! input, eigen energies and occupations of noninteracting part
  real(dp), dimension(norb, norb) :: Skq
  ! input, Structural function, Skq(ibnd, jbnd)
  !   ibnd is band index for k; jbnd is band index for kq
  !
  integer ibnd, jbnd
  real(dp) :: occ_diff
  !
  chi=cmplx_0
  !
  do ibnd=1, norb
    do jbnd=1, norb
      !
      occ_diff=occ_kq(jbnd)-occ_k(ibnd)
      if (abs(occ_diff)>eps9) then
        !
        chi=chi-Skq(ibnd, jbnd)*occ_diff/(w+eig_kq(jbnd)-eig_k(ibnd))
        !
      endif
      !
    enddo ! jbnd
  enddo ! ibnd
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_trace(chi0, nu, nnu, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   of frequency nu
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, beta_in
  use para,         only : inode, para_merge_cmplx, distribute_calc, first_idx, last_idx
  !
  implicit none
  !
  integer nnu
  complex(dp), dimension(nnu) :: chi0, nu
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  complex(dp), dimension(norb, norb) :: hk, hkq
  real(dp), dimension(norb) :: ek, ekq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
  !
  complex(dp), dimension(nnu) :: chi0_part
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
    call calc_eig_occ_Skq(ek, ekq, occ_k, occ_kq, Skq, hk, hkq, beta_in)
    !
    do inu=1, nnu
      call calc_partial_chi_bare_trace(chi0_part(inu), nu(inu), ek, ekq, occ_k, occ_kq, Skq)
    enddo
    chi0=chi0+chi0_part
    !
  enddo
  !
  call para_merge_cmplx(chi0, nnu)
  !
  chi0(:)=chi0(:)/nkpt
  !
END SUBROUTINE

SUBROUTINE calc_chi_corr_trace(chi, chi0, nnu, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   Up to first nnu Matsubara frequencies
  !
  use constants,    only : dp, cmplx_0, stdout, eps6
  use wanndata,     only : norb
  use input,        only : nkx, nky, nkz, nexact
  use impurity,     only : beta
  use para,         only : inode
  !
  implicit none
  !
  integer      :: nnu
  complex(dp), dimension(nnu) :: chi, chi0
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  complex(dp), dimension(norb, norb) :: hk, hkq
  real(dp), dimension(norb) :: ek, ekq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
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
        if (inode.eq.0) then
          write(stdout, '(A,3F9.4)') ' Working on k: ', k
        endif
        !
        call calc_hk(hk, k)
        call calc_hk(hkq, kq)
        !
        call calc_eig_occ_Skq(ek, ekq, occ_k, occ_kq, Skq, hk, hkq, beta)
        !
        call calc_partial_chi_corr_trace(chi_part, chi0_part, nnu, nexact, hk, hkq, ek, ekq, occ_k, occ_kq, Skq)
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

SUBROUTINE calc_partial_chi_corr_trace(chi, chi0, nnu, nexact, hk, hkq, eig_k, eig_kq, occ_k, occ_kq, Skq)
  !
  ! This subroutine calculates part of the susceptibility (one kpt in the K-mesh),
  !   Parallel over frequency
  !
  use constants, only : dp, eps9, cmplx_0, twopi, cmplx_i, stdout
  use wanndata,  only : norb
  use impurity,  only : beta
  use para,      only : distribute_calc, first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx
  !
  implicit none
  !
  integer       :: nnu
  ! First nnu Matsubara frequencies will be calculated
  integer       :: nexact
  ! We only calculate -nexact to nexact contributions exactly (the difference decays very fast)
  complex(dp), dimension(nnu) :: chi, chi0
  ! chi: correlated susceptibility
  ! chi0: noninteracting susceptibility
  !
  complex(dp), dimension(norb, norb) :: hk, hkq
  ! non-interacting Hamiltonian at k and k+q
  real(dp), dimension(norb) :: eig_k, eig_kq, occ_k, occ_kq
  ! non-interacting eigen energies and occupations
  real(dp), dimension(norb, norb) :: Skq
  ! structural function in non-interacting case
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
    call calc_partial_chi_bare_trace(chi0(inu+1), nu, eig_k, eig_kq, occ_k, occ_kq, Skq)
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
          chi_diff=gk(ibnd, jbnd, iom-first_idx+1)*gkq(jbnd, ibnd, iom-first_idx+1)-Skq(ibnd, jbnd)/((w1-eig_k(ibnd))*(w2-eig_kq(jbnd)))
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

SUBROUTINE analyze_part_chi_corr_trace(chi, chi0, nnu, hk, hkq, eig_k, eig_kq, occ_k, occ_kq, Skq)
  !
  !  This subroutine analyzes the composition of chi_corr
  !    Please be noted that unlike the calculation subroutine
  !      nnu here is the matsubara frequency being analyzed
  !      and chi/chi0 contains the frequency components instead of final results
  !
  use constants, only : dp, cmplx_0, twopi, cmplx_i
  use wanndata,  only : norb
  use impurity,  only : nfreq, beta
  use para,      only : first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx, distribute_calc
  !
  implicit none
  !
  integer       :: nnu
  complex(dp), dimension(2*nfreq-nnu) :: chi, chi0
  complex(dp), dimension(norb)  :: hk, hkq
  real(dp), dimension(norb)     :: eig_k, eig_kq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
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
        chi0(iom)=chi0(iom)+Skq(ibnd, jbnd)/((w1-eig_k(ibnd))*(w2-eig_kq(jbnd)))
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

SUBROUTINE analyze_part_chi_corr_trace1(chi, chif, chic, nnu, hk, hkq, eig_k, eig_kq, occ_k, occ_kq, Skq)
  !
  !  This subroutine analyzes the composition of chi_corr
  !    Please be noted that unlike the calculation subroutine
  !      nnu here is the matsubara frequency being analyzed
  !      and chi/chi0 contains the frequency components instead of final results
  !
  use constants, only : dp, cmplx_0, twopi, cmplx_i
  use wanndata,  only : norb
  use impurity,  only : nfreq, beta, fulldim, basis_map
  use para,      only : first_idx, last_idx, para_merge_cmplx, inode, para_collect_cmplx, para_distribute_cmplx, distribute_calc
  !
  implicit none
  !
  integer       :: nnu
  complex(dp), dimension(2*nfreq-nnu) :: chi, chif, chic
  complex(dp), dimension(norb)  :: hk, hkq
  real(dp), dimension(norb)     :: eig_k, eig_kq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
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
  use input,        only : nkx, nky, nkz
  use impurity,     only : beta, nfreq
  use para,         only : inode
  !
  implicit none
  !
  integer     :: nnu
  !
  complex(dp), dimension(2*nfreq-nnu) :: chi, chi0
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  complex(dp), dimension(norb, norb) :: hk, hkq
  real(dp), dimension(norb) :: ek, ekq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
  !
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
        call calc_eig_occ_Skq(ek, ekq, occ_k, occ_kq, Skq, hk, hkq, beta)
        !
        call analyze_part_chi_corr_trace(chi_part, chi0_part, nnu, hk, hkq, ek, ekq, occ_k, occ_kq, Skq)
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
  use input,        only : nkx, nky, nkz
  use impurity,     only : beta, nfreq
  use para,         only : inode
  !
  implicit none
  !
  integer     :: nnu
  !
  complex(dp), dimension(2*nfreq-nnu) :: chi, chif, chic
  real(dp), dimension(3) :: qv
  !
  real(dp), dimension(3) :: k, kq
  complex(dp), dimension(norb, norb) :: hk, hkq
  real(dp), dimension(norb) :: ek, ekq, occ_k, occ_kq
  real(dp), dimension(norb, norb) :: Skq
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
        call calc_eig_occ_Skq(ek, ekq, occ_k, occ_kq, Skq, hk, hkq, beta)
        !
        call analyze_part_chi_corr_trace1(chi_part, chif_part, chic_part, nnu, hk, hkq, ek, ekq, occ_k, occ_kq, Skq)
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
