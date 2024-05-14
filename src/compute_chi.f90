MODULE chi_internal
  !
  use constants,  only : dp
  !
  implicit none
  !
  ! THESE VARIABLES ARE ALWAYS USED IF FAST CALCULATION IS REQUESTED
  !
  integer, dimension(:),  allocatable :: kq_idx
  ! THESE VARIABLES ARE USED IF CHI0 IS CALCULATED USING LEHMAN REP
  !
  real(dp), dimension(:), allocatable :: ek, ekq, occ_k, occ_kq
  ! Eigen value and occupation at k and k+q
  !
  complex(dp), dimension(:, :), allocatable :: hk, hkq
  ! Hamiltonian / Eigen state at k and k+q
  !
  real(dp), dimension(:, :), allocatable :: Skq
  ! Form factor between k & k+q
  !   Skq(ibnd, jbnd)
  !   ibnd is band index for k; jbnd is band index for kq
  !   Facilitates fast trace
  !
  ! THESE VARIABLES ARE USED IF LEHMAN FAST CALCULATION IS REQUESTED
  !
  real(dp), dimension(:, :), allocatable       :: eig_full, occ_full
  complex(dp), dimension(:, :, :), allocatable :: Uk_full
  complex(dp), dimension(:, :, :), allocatable :: Uk_local, Ukq_local
  !
  ! THESE VARIABLES ARE USED IF G*G ALGORITHM IS EMPLOYED (for correlated case)
  !
  complex(dp), dimension(:, :), allocatable :: Hff_k, Hff_kq
  complex(dp), dimension(:, :), allocatable :: Vfc_k, Vfc_kq
  real(dp), dimension(:), allocatable       :: Ecc_k, Ecc_kq
  !
  ! THESE VARIABLES ARE USED IF G*G FAST ALGORITHM IS EMPLOYED (for correlated case)
  !
  complex(dp), dimension(:, :, :), allocatable :: Hff_full, Vfc_full, Hkff_loc, Hkqff_loc, Vkfc_loc, Vkqfc_loc
  real(dp), dimension(:, :), allocatable       :: Ecc_full, Ecck_loc, Ecckq_loc
  !
  ! THESE VARIABLES ARE USED IN CORRELATED CASE
  !
  complex(dp), dimension(:, :, :), allocatable :: Gff_k, Gff_kq
  !
  ! THESE VARIABLES ARE USED FOR SUPERFAST ALGORITHM IN CORRELATED CASE
  !
  complex(dp), dimension(:, :, :, :), allocatable :: Gff_full, Gkff_loc, Gkqff_loc
  !
  !
  ! THESE VARIABLES ARE USED FOR FULL MATRIX CALCULATIONS
  !
  complex(dp), dimension(:, :, :), allocatable :: chiff, chicc, chifc, chicf
  !
  ! THESE VARIABLES ARE USED FOR RPA CALCULATIONS
  complex(dp), dimension(:, :, :), allocatable :: chi0ff, chi0cc, chi0fc, chi0cf
  !
  contains
  !
  subroutine finalize_chi_internal()
    !
    use constants,    only : fout3, fout4, fout5, fout6, fin3, fin4, fin5, fin6
    !
    implicit none
    !
    logical isopen
    !
    if (allocated(chiff))         deallocate(chiff)
    if (allocated(chicc))         deallocate(chicc)
    if (allocated(chifc))         deallocate(chifc)
    if (allocated(chicf))         deallocate(chicf)
    if (allocated(chi0ff))        deallocate(chi0ff)
    if (allocated(chi0cc))        deallocate(chi0cc)
    if (allocated(chi0fc))        deallocate(chi0fc)
    if (allocated(chi0cf))        deallocate(chi0cf)
    !
    if (allocated(ek))            deallocate(ek)
    if (allocated(ekq))           deallocate(ekq)
    if (allocated(occ_k))         deallocate(occ_k)
    if (allocated(occ_kq))        deallocate(occ_kq)
    if (allocated(hk))            deallocate(hk)
    if (allocated(hkq))           deallocate(hkq)
    if (allocated(Skq))           deallocate(Skq)
    !
    if (allocated(kq_idx))        deallocate(kq_idx)
    if (allocated(eig_full))      deallocate(eig_full)
    if (allocated(occ_full))      deallocate(occ_full)
    if (allocated(Uk_full))       deallocate(Uk_full)
    if (allocated(Uk_local))      deallocate(Uk_local)
    if (allocated(Ukq_local))     deallocate(Ukq_local)
    !
    if (allocated(Hff_full))      deallocate(Hff_full)
    if (allocated(Hff_k))         deallocate(Hff_k)
    if (allocated(hff_kq))        deallocate(Hff_kq)
    if (allocated(Vfc_full))      deallocate(Vfc_full)
    if (allocated(Vfc_k))         deallocate(Vfc_k)
    if (allocated(Vfc_kq))        deallocate(Vfc_kq)
    if (allocated(Ecc_full))      deallocate(Ecc_fulL)
    !
    inquire(unit=fin3, opened=isopen)
    if (isopen) close(unit=fin3)
    !
    inquire(unit=fin4, opened=isopen)
    if (isopen) close(unit=fin4)
    !
    inquire(unit=fin5, opened=isopen)
    if (isopen) close(unit=fin5)
    !
    inquire(unit=fin6, opened=isopen)
    if (isopen) close(unit=fin6)
    !
    inquire(unit=fout3, opened=isopen)
    if (isopen) close(unit=fout3)
    !
    inquire(unit=fout4, opened=isopen)
    if (isopen) close(unit=fout4)
    !
    inquire(unit=fout5, opened=isopen)
    if (isopen) close(unit=fout5)
    !
    inquire(unit=fout6, opened=isopen)
    if (isopen) close(unit=fout6)
    !
  end subroutine

  subroutine init_chi_matrix(nw)
    !
    use constants, only : dp, fout3, fout4, fout5, fout6
    use lattice,   only : ham, nkirr
    use para,      only : first_idx, last_idx, inode
    use IntRPA,    only : nFFidx, nCCidx
    !
    implicit none
    !
    integer nw
    !
    allocate(chiff(nFFidx, nFFidx, nw))
    !
    if (nCCidx>0) then
      !
      allocate(chifc(nFFidx, nCCidx, nw))
      allocate(chicf(nCCidx, nFFidx, nw))
      allocate(chicc(nCCidx, nCCidx, nw))
      !
    endif
    !
    if (inode.eq.0) then
      !
      open(unit=fout3, file='chiff.dat', form='unformatted', access='direct', recl=nFFidx*nFFidx*16)
      if (nCCidx>0) then
        open(unit=fout4, file='chicc.dat', form='unformatted', access='direct', recl=nCCidx*nCCidx*16)
        open(unit=fout5, file='chifc.dat', form='unformatted', access='direct', recl=nFFidx*nCCidx*16)
        open(unit=fout6, file='chicf.dat', form='unformatted', access='direct', recl=nCCidx*nFFidx*16)
      endif
      !
    endif
    !
  end subroutine

  subroutine init_chi_matrix_RPA(nw)
    !
    use constants, only : dp, fin3, fin4, fin5, fin6, fout3, fout4, fout5, fout6
    use lattice,   only : ham, nkirr
    use para,      only : first_idx, last_idx, inode
    use IntRPA,    only : nFFidx, nCCidx
    !
    implicit none
    !
    integer nw
    !
    allocate(chiff(nFFidx, nFFidx, nw))
    allocate(chi0ff(nFFidx, nFFidx, nw))
    !
    if (nCCidx>0) then
      !
      allocate(chifc(nFFidx, nCCidx, nw))
      allocate(chicf(nCCidx, nFFidx, nw))
      allocate(chicc(nCCidx, nCCidx, nw))
      allocate(chi0fc(nFFidx, nCCidx, nw))
      allocate(chi0cf(nCCidx, nFFidx, nw))
      allocate(chi0cc(nCCidx, nCCidx, nw))
      !
    endif
    !
    if (inode.eq.0) then
      !
      open(unit=fout3, file='chiff.dat', form='unformatted', access='direct', recl=nFFidx*nFFidx*16)
      open(unit=fin3, file='chi0ff.dat', form='unformatted', access='direct', recl=nFFidx*nFFidx*16)
      if (nCCidx>0) then
        open(unit=fout4, file='chicc.dat', form='unformatted', access='direct', recl=nCCidx*nCCidx*16)
        open(unit=fout5, file='chifc.dat', form='unformatted', access='direct', recl=nFFidx*nCCidx*16)
        open(unit=fout6, file='chicf.dat', form='unformatted', access='direct', recl=nCCidx*nFFidx*16)
        open(unit=fin4, file='chi0cc.dat', form='unformatted', access='direct', recl=nCCidx*nCCidx*16)
        open(unit=fin5, file='chi0fc.dat', form='unformatted', access='direct', recl=nFFidx*nCCidx*16)
        open(unit=fin6, file='chi0cf.dat', form='unformatted', access='direct', recl=nCCidx*nFFidx*16)
      endif
      !
    endif
    !
  end subroutine

  subroutine init_chi_internal_lehman(fast_calc, step)
    !
    use constants, only : dp
    use lattice,   only : ham, nkirr
    use para,      only : first_idx, last_idx, inode
    !
    implicit none
    !
    integer step
    logical fast_calc
    !
    ! These are always required in any cases...
    !
    if (step.ne.2) then
      allocate(ek(ham%norb), ekq(ham%norb), occ_k(ham%norb), occ_kq(ham%norb))
      allocate(hk(ham%norb, ham%norb), hkq(ham%norb, ham%norb), Skq(ham%norb, ham%norb))
    endif
    !
    if (fast_calc) then
      !
      if (step==1) then
        !
        allocate(eig_full(ham%norb, nkirr), occ_full(ham%norb, nkirr))
        allocate(kq_idx(nkirr))
        allocate(Uk_local(ham%norb, ham%norb, last_idx-first_idx+1), Ukq_local(ham%norb, ham%norb, last_idx-first_idx+1))
        !
      else
        !
        if (inode.eq.0) then
          allocate(Uk_full(ham%norb, ham%norb, nkirr))
        else
          allocate(Uk_full(1, 1, 1))
        endif
        !
      endif
      !
    endif
    !
  end subroutine

  subroutine init_chi_internal_GG(fast_calc, step)
    !
    use constants, only : dp
    use lattice,   only : ham, nkirr
    use para,      only : first_idx, last_idx, inode
    !
    implicit none
    !
    logical fast_calc
    integer step
    !
    ! These are always required in any cases...
    !
    if (step.ne.2) then
      allocate(ek(ham%norb), ekq(ham%norb))
      allocate(hk(ham%norb, ham%norb), hkq(ham%norb, ham%norb))
    endif
    !
    if (fast_calc) then
      !
      if (step==1) then
        !
        allocate(kq_idx(nkirr))
        allocate(eig_full(ham%norb, nkirr))
        allocate(Uk_local(ham%norb, ham%norb, last_idx-first_idx+1), Ukq_local(ham%norb, ham%norb, last_idx-first_idx+1))
        !
      else
        !
        if (inode.eq.0) then
          allocate(Uk_full(ham%norb, ham%norb, nkirr))
        else
          allocate(Uk_full(1, 1, 1))
        endif
        !
      endif
      !
    endif
    !
  end subroutine

  subroutine init_chi_internal_corr(fast_calc, step)
    !
    use constants,    only : dp
    use lattice,      only : ham, ndimf, ndimc, nkirr
    use para,         only : first_idx, last_idx, inode
    use input,        only : nnu, npade
    use IntRPA,       only : nFFidx, nCCidx
    !
    implicit none
    !
    logical fast_calc
    integer step
    !
    ! Always needed
    !
    if (step.ne.2) then
      !
      allocate(Gff_k(ndimf, ndimf, 2*npade), Gff_kq(ndimf, ndimf, 2*npade))
      !
    endif
    !
    if (.not. fast_calc) then
      !
      if (step==1) then
        !
        allocate(kq_idx(nkirr))
        !
        allocate(Hff_k(ndimf, ndimf), Vfc_k(ndimf, ndimc), Ecc_k(ndimc))
        allocate(Hff_kq(ndimf, ndimf), Vfc_kq(ndimf, ndimc), Ecc_kq(ndimc))
        !
        allocate(Ecc_full(ndimc,        nkirr))
        allocate(Hkff_loc(ndimf, ndimf, last_idx-first_idx+1))
        allocate(Vkfc_loc(ndimf, ndimc, last_idx-first_idx+1))
        allocate(Ecck_loc(ndimc,        last_idx-first_idx+1))
        !
        allocate(Hkqff_loc(ndimf, ndimf, last_idx-first_idx+1))
        allocate(Vkqfc_loc(ndimf, ndimc, last_idx-first_idx+1))
        allocate(Ecckq_loc(ndimc,        last_idx-first_idx+1))
        !
      else
        !
        if (inode.eq.0) then
          !
          allocate(Hff_full(ndimf, ndimf, nkirr))
          allocate(Vfc_full(ndimf, ndimc, nkirr))
          !
        else
          !
          allocate(Hff_full(1, 1, 1))
          allocate(Vfc_full(1, 1, 1))
          !
        endif
        !
      endif
      !
    else
      !
      if (step==1) then
        !
        allocate(kq_idx(nkirr))
        !
        allocate(Hkff_loc(ndimf, ndimf, last_idx-first_idx+1))
        allocate(Vkfc_loc(ndimf, ndimc, last_idx-first_idx+1))
        allocate(Ecck_loc(ndimc,        last_idx-first_idx+1))
        !
      else
        !
        allocate(Gkff_loc(ndimf, ndimf, 2*npade, last_idx-first_idx+1))
        allocate(Gkqff_loc(ndimf, ndimf, 2*npade, last_idx-first_idx+1))
        !
        if (inode.eq.0) then
          !
          allocate(Gff_full(ndimf, ndimf, 2*npade, nkirr))
          !
        else
          !
          allocate(Gff_full(1, 1, 1, 1))
          !
        endif
        !
      endif
      !
    endif
    !
  end subroutine

  subroutine save_chi_matrix_one(idx, ii, ff_only)
    !
    use constants, only : fout3, fout4, fout5, fout6
    use para,      only : inode
    use IntRPA,    only : nCCidx, nFFidx
    !
    implicit none
    !
    integer idx
    logical ff_only
    integer ii
    !
    if (inode.eq.0) then
      write(fout3, rec=idx+ii) chiff(:, :, ii)
      !
      if ((.not. ff_only).and.(nCCidx.ne.0)) then
        write(fout4, rec=idx+ii) chicc(:, :, ii)
        write(fout5, rec=idx+ii) chifc(:, :, ii)
        write(fout6, rec=idx+ii) chicf(:, :, ii)
      endif
      !
    endif
    !
  end subroutine

  subroutine save_chi_matrix(idx, nn, ff_only)
    !
    use constants, only : fout3, fout4, fout5, fout6
    use para,      only : inode
    use IntRPA,    only : nCCidx, nFFidx
    !
    implicit none
    !
    integer idx, nn
    logical ff_only
    integer ii
    !
    if (inode.eq.0) then
      do ii=1, nn
        write(fout3, rec=idx+ii) chiff(:, :, ii)
        !
        if ((.not. ff_only).and.(nCCidx.ne.0)) then
          write(fout4, rec=idx+ii) chicc(:, :, ii)
          write(fout5, rec=idx+ii) chifc(:, :, ii)
          write(fout6, rec=idx+ii) chicf(:, :, ii)
        endif
        !
      enddo
    endif
    !
  end subroutine

  subroutine read_chi_matrix(idx, nn, ff_only)
    !
    use constants, only : fout3, fout4, fout5, fout6
    use para,      only : inode
    use IntRPA,    only : nCCidx, nFFidx
    !
    implicit none
    !
    integer idx, nn
    logical ff_only
    integer ii
    !
    if (inode.eq.0) then
      do ii=1, nn
        read(fout3, rec=idx+ii) chiff(:, :, ii)
        !
        if ((.not. ff_only).and.(nCCidx.ne.0)) then
          read(fout4, rec=idx+ii) chicc(:, :, ii)
          read(fout5, rec=idx+ii) chifc(:, :, ii)
          read(fout6, rec=idx+ii) chicf(:, :, ii)
        endif
        !
      enddo
    endif
    !
  end subroutine

  subroutine read_chi_matrix_RPA(idx, nn, ff_only)
    !
    use constants, only : fin3, fin4, fin5, fin6
    use para,      only : inode
    use IntRPA,    only : nCCidx, nFFidx
    !
    implicit none
    !
    integer idx, nn
    logical ff_only
    integer ii
    !
    if (inode.eq.0) then
      do ii=1, nn
        read(fin3, rec=idx+ii) chi0ff(:, :, ii)
        !
        if ((.not. ff_only).and.(nCCidx.ne.0)) then
          read(fin4, rec=idx+ii) chi0cc(:, :, ii)
          read(fin5, rec=idx+ii) chi0fc(:, :, ii)
          read(fin6, rec=idx+ii) chi0cf(:, :, ii)
        endif
        !
      enddo
    endif
    !
  end subroutine

  subroutine show_chi_diag(iou, iw)
    !
    use IntRPA,    only : nFFidx, nCCidx
    use input,     only : ff_only
    !
    implicit none
    !
    integer iou, iw
    integer ii
    !
    if (.not. allocated(chiff)) then
      !
      write(*, *) "!!! INTERNAL ERROR, requires full matrix?"
      stop
      !
    endif
    !
    write(iou, *) "  Diagonal elements in FF:"
    write(iou, *) "    Real Part:"
    !
    do ii=1, nFFidx
      write(iou, '(1F14.9)', advance='no') real(chiff(ii, ii, iw))
      if (MOD(ii,10)==0.or.ii==nFFidx) write(iou, *) 
    enddo
    write(iou, *)
    write(iou, *) "    Imag Part:"
    do ii=1, nFFidx
      write(iou, '(1F14.9)', advance='no') aimag(chiff(ii, ii, iw))
      if (MOD(ii,10)==0.or.ii==nFFidx) write(iou, *) 
    enddo
    write(iou, *)
    !
    if (nCCidx>0.and.(.not.ff_only)) then
      !
      write(iou, *) "  Diagonal elements in CC:"
      write(iou, *) "    Real Part:"
      do ii=1, nCCidx
        write(iou, '(1F14.9)', advance='no') real(chicc(ii, ii, iw))
        if (MOD(ii,10)==0.or.ii==nCCidx) write(iou, *) 
      enddo
      write(iou, *)
      write(iou, *) "    Imag Part:"
      do ii=1, nCCidx
        write(iou, '(1F14.9)', advance='no') aimag(chicc(ii, ii, iw))
        if (MOD(ii,10)==0.or.ii==nCCidx) write(iou, *) 
      enddo
      write(iou, *)
      !
    endif
    !
  end subroutine
  !
END MODULE

SUBROUTINE prepare_lehman
  !
  use constants,  only : dp, stdout
  use input,      only : beta, nnu, fast_calc, trace_only
  use lattice,    only : nkirr, ham, kvec
  use wanndata,   only : calc_hk, finalize_wann
  use para,       only : distribute_calc, first_idx, last_idx, para_collect_cmplx, para_merge_real, inode, nnode
  use linalgwrap, only : eigen
  use chi_internal, only : hk, eig_full, occ_full, Uk_full, Uk_local, init_chi_internal_lehman, init_chi_matrix
  use IntRPA,     only : nFFidx, nCCidx
  !
  implicit none
  !
  integer ik, ii
  real(dp)      :: mem_req
  !
  real(dp) calc_occ
  !
  if (.not. trace_only) then
    !
    call init_chi_matrix(nnu)
    !
  endif
  !
  if (inode.eq.0) then
    !
    mem_req = (ham%norb * ham%norb * ham%nrpt * nnode * 16.d0)
    !
    write(stdout, '(A)') ' Lehman representation calculation '
    write(stdout, '(A)') '   This is faster than G*G algorithm, '
    write(stdout, '(A)') '   but may NOT be accurate for highly degenerate systems. '
    write(stdout, '(A)') 
    write(stdout, '(A,1F9.1,A)') ' Wannier Hamiltonian takes ~ ', mem_req/1.d9, 'GB'
    !
    if (fast_calc) then
      !
      if (nkirr > ham%nrpt * nnode) then
        mem_req = (ham%norb * ham%norb * nkirr * 2.d0) *16.d0
      else
        mem_req = mem_req + (ham%norb * ham%norb * nkirr) *16.d0
      endif
      !
      write(stdout, *) '  !!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!'
      write(stdout, *) '  FAST ALGORITHM, MAKE SURE YOU HAVE ENOUGH MEMORY'
      write(stdout, *) '    AND qv IS COMMENSURATE!!!'
      !
      if (.not. trace_only) then
        mem_req = mem_req + (nFFidx + nCCidx) * (nFFidx + nCCidx) * nnu * nnode * 16.d0
      endif
      !
      write(stdout, '(A,1F9.1,A)') '    Estimated required memory size:', mem_req/1.d9, 'GB'
      write(stdout, *) '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      !
    endif
    !
  endif
  !
  if (fast_calc) then
    ! Do Full Kmesh diagonalization beforehand
    !    in case of fast calculation
    !
    call distribute_calc(nkirr)
    ! Always parallel in Kmesh in fast calculations
    !
    call init_chi_internal_lehman(fast_calc, 1)
    !
    eig_full=0.0
    occ_full=0.0
    !
    do ik=first_idx, last_idx
      !
      call calc_hk(hk, ham, kvec(:, ik))
      call eigen(eig_full(:, ik), hk, ham%norb)
      !
      Uk_local(:, :, ik-first_idx+1)=hk(:, :)
      !
      do ii=1, ham%norb
        occ_full(ii, ik)=calc_occ(eig_full(ii, ik), beta)
      enddo
      !
    enddo
    !
    call finalize_wann(ham, .false.)
    call init_chi_internal_lehman(fast_calc, 2)
    !
    call para_merge_real(eig_full, ham%norb*nkirr)
    call para_merge_real(occ_full, ham%norb*nkirr)
    call para_collect_cmplx(Uk_full, Uk_local, ham%norb*ham%norb)
    !
  else
    !
    call init_chi_internal_lehman(fast_calc, 1)
    !
  endif
  !
END SUBROUTINE

SUBROUTINE prepare_GG
  !
  use constants,  only : dp, stdout
  use input,      only : beta, nnu, fast_calc, trace_only, npade
  use lattice,    only : nkirr, ham, kvec
  use wanndata,   only : calc_hk, finalize_wann
  use para,       only : distribute_calc, first_idx, last_idx, para_collect_cmplx, para_merge_real, inode, nnode
  use linalgwrap, only : eigen
  use chi_internal, only : hk, eig_full, Uk_full, Uk_local, init_chi_internal_GG, init_chi_matrix
  use pade_sum,   only : init_pade
  use IntRPA,     only : nCCidx, nFFidx
  !
  implicit none
  !
  integer ik, ii
  real(dp) mem_req
  !
  call init_chi_matrix(nnu)
  !
  call init_pade(npade, .true.)
  !
  if (inode.eq.0) then
    !
    mem_req = (ham%norb * ham%norb * ham%nrpt * nnode * 16.d0)
    !
    write(stdout, '(A,1F9.1,A)') ' Wannier Hamiltonian takes ~ ', mem_req/1.d9, 'GB'
    !
    if (fast_calc) then
      !
      if (nkirr > ham%nrpt * nnode) then
        mem_req = (ham%norb * ham%norb * nkirr * 2.d0) *16.d0
      else
        mem_req = mem_req + (ham%norb * ham%norb * nkirr) *16.d0
      endif
      !
      write(stdout, *) '  !!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!'
      write(stdout, *) '  FAST ALGORITHM, MAKE SURE YOU HAVE ENOUGH MEMORY'
      write(stdout, *) '    AND qv IS COMMENSURATE!!!'
      !
      mem_req = mem_req + (nFFidx + nCCidx) * (nFFidx + nCCidx) * nnu * nnode * 16.d0
      !
      write(stdout, '(A,1F9.1,A)') '    Estimated required memory size:', mem_req/1.d9, 'GB'
      write(stdout, *) '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      !
    endif
    !
  endif
  !
  if (fast_calc) then
    ! Do Full Kmesh diagonalization beforehand
    !    in case of fast calculation
    !
    call distribute_calc(nkirr)
    ! Always parallel in Kmesh in fast calculations
    !
    call init_chi_internal_GG(fast_calc, 1)
    !
    eig_full=0.0
    !
    do ik=first_idx, last_idx
      !
      call calc_hk(hk, ham, kvec(:, ik))
      call eigen(eig_full(:, ik), hk, ham%norb)
      !
      Uk_local(:, :, ik-first_idx+1)=hk(:, :)
      !
    enddo
    !
    call finalize_wann(ham, .false.)
    call init_chi_internal_GG(fast_calc, 2)
    !
    call para_merge_real(eig_full, ham%norb*nkirr)
    call para_collect_cmplx(Uk_full, Uk_local, ham%norb*ham%norb)
    !
  else
    !
    call init_chi_internal_GG(fast_calc, 1)
    !
  endif
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_trace_lehman_kernel(chi0, w, nw)
  !
  ! This subroutine calculates part of the susceptibility,
  !   a summation over k (integration over BZ) is required to obtain full chi
  !
  use constants, only : dp, eps12, cmplx_0
  use chi_internal, only : ek, ekq, Skq, occ_k, occ_kq
  use lattice,   only : ham
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  ! chi0: output susceptibility
  !   w: input frequency
  !
  integer ibnd, jbnd
  real(dp) :: occ_diff
  !
  do ibnd=1, ham%norb
    do jbnd=1, ham%norb
      !
      occ_diff=occ_k(ibnd)-occ_kq(jbnd)
      if (abs(occ_diff)>eps12) then
        !
        chi0(:)=chi0(:)+Skq(ibnd, jbnd)*occ_diff/(w(:)+ekq(jbnd)-ek(ibnd))
        !
      endif
      !
    enddo ! jbnd
  enddo ! ibnd
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_trace_lehman(chi0, w, nw, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   of frequency w
  !
  use constants,    only : dp, cmplx_0, stdout, fdebug
  use input,        only : beta
  use wanndata,     only : calc_hk
  use lattice,      only : nkirr, kvec, ham
  use chi_internal, only : ek, ekq, occ_k, occ_kq, hk, hkq, Skq
  use para,         only : para_merge_cmplx, distribute_calc, first_idx, last_idx!, inode
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  real(dp),    dimension(3)  :: kq
  integer                    :: ik, iw
  !
  !character(len=40) :: fn
  !
  chi0=cmplx_0
  !
  call distribute_calc(nkirr)
  !
  !if (inode<10) then
  !  write(fn, "(A9,I1)"), "DEBUGOUT_", inode
  !elseif (inode<100) then
  !  write(fn, "(A9,I2)"), "DEBUGOUT_", inode
  !elseif (inode<1000) then
  !  write(fn, "(A9,I3)"), "DEBUGOUT_", inode
  !endif

  !open(unit=fdebug+inode, file=trim(fn))
  !
  do ik=first_idx, last_idx
    !
    !write(fdebug+inode, *) "ik: ", ik
    !write(fdebug+inode, *) "kvec: ", kvec(:, ik)
    call calc_hk(hk, ham, kvec(:, ik))
    call calc_eig_occ_U(ek, occ_k, hk, beta)
    !
    !write(fdebug+inode, *) "ek, occ:"
    !write(fdebug+inode, '(5F9.4)') ek
    !write(fdebug+inode, '(5F9.4)') occ_k
    !write(fdebug+inode, *) "egv:"
    !call show_matrix(hk, ham%norb, fdebug+inode)
    !
    kq=kvec(:, ik)+qv(:)
    call calc_hk(hkq, ham, kq)
    call calc_eig_occ_U(ekq, occ_kq, hkq, beta)
    !
    !write(fdebug+inode, *) "ek, occ:"
    !write(fdebug+inode, '(5F9.4)') ekq
    !write(fdebug+inode, '(5F9.4)') occ_kq
    !write(fdebug+inode, *) "egv:"
    !call show_matrix(hkq, ham%norb, fdebug+inode)
    !
    call calc_Skq
    !
    !write(fdebug+inode, *) " Skq: "
    !call show_rmatrix(Skq, ham%norb, fdebug+inode)
    !
    call calc_chi_bare_trace_lehman_kernel(chi0, w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chi0, nw)
  !
  chi0(:)=chi0(:)/nkirr
  !
  !close(fdebug+inode)
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_trace_lehman_fast(chi0, w, nw, qv)
  !
  ! This subroutine calculates susceptibility at qv
  !   of frequency w
  !
  use constants,    only : dp, cmplx_0, twopi, fdebug
  use lattice,      only : nkirr, kvec, ham
  use para,         only : para_merge_cmplx, first_idx, last_idx !, inode
  use chi_internal, only : ek, ekq, Skq, occ_k, occ_kq, eig_full, occ_full, hk, hkq, Uk_local, Ukq_local, kq_idx
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  integer                    :: ik, ii, jj
  real(dp), dimension(3, last_idx-first_idx+1)  :: kq_local
  real(dp),    dimension(3)  :: DG
  real(dp)                   :: GdotTau
  complex(dp), dimension(ham%norb) :: extra_phase
  complex(dp)                :: fact
  !character(len=40) :: fn
  !
  chi0=cmplx_0

  !if (inode<10) then
  !  write(fn, "(A9,I1)"), "DEBUGOUT_", inode
  !elseif (inode<100) then
  !  write(fn, "(A9,I2)"), "DEBUGOUT_", inode
  !elseif (inode<1000) then
  !  write(fn, "(A9,I3)"), "DEBUGOUT_", inode
  !endif

  !open(unit=fdebug+inode, file=trim(fn))
  !
  do ik=first_idx, last_idx
    !
    kq_local(:, ik-first_idx+1)=kvec(:, ik)+qv(:)
    !
  enddo
  ! Locate k+q elements in node 0 and populate
  call find_idx_in_kmesh(kq_local)
  call reorder_and_distribute_Ukq()
  !
  do ik=first_idx, last_idx
    !
    !write(fdebug+inode, *) "ik: ", ik
    !write(fdebug+inode, *) "kvec: ", kvec(:, ik)
    ek(:)=eig_full(:, ik)
    ekq(:)=eig_full(:, kq_idx(ik))
    !
    occ_k(:)=occ_full(:, ik)
    occ_kq(:)=occ_full(:, kq_idx(ik))
    !
    hk(:, :)=Uk_local(:, :, ik-first_idx+1)
    hkq(:, :)=Ukq_local(:, :, ik-first_idx+1)
    !
    !write(fdebug+inode, *) "ek, occ:"
    !write(fdebug+inode, '(5F9.4)') ek
    !write(fdebug+inode, '(5F9.4)') occ_k
    !write(fdebug+inode, *) "egv:"
    !call show_matrix(hk, ham%norb, fdebug+inode)
    !
    !write(fdebug+inode, *) "ek, occ:"
    !write(fdebug+inode, '(5F9.4)') ekq
    !write(fdebug+inode, '(5F9.4)') occ_kq
    !write(fdebug+inode, *) "egv:"
    !call show_matrix(hkq, ham%norb, fdebug+inode)

    ! Fix the phase due to site position
    DG(:)=kq_local(:, ik-first_idx+1)-kvec(:, kq_idx(ik))
    !
    do ii=1, ham%norb
      !
      GdotTau=SUM(DG(:)*ham%tau(:, ii))*twopi
      extra_phase(ii)=CMPLX(cos(GdotTau), sin(GdotTau), KIND=dp)
      !
    enddo
    !
    do ii=1, ham%norb
      do jj=1, ham%norb
        fact=sum( extra_phase(:)*hk(:, ii)*conjg(hkq(:, jj)) )
        Skq(ii, jj)=conjg(fact)*fact
      enddo
    enddo
    !
    !write(fdebug+inode, *) " Skq: "
    !call show_rmatrix(Skq, ham%norb, fdebug+inode)
    !
    call calc_chi_bare_trace_lehman_kernel(chi0, w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chi0, nw)
  !
  chi0(:)=chi0(:)/nkirr
  !
  !close(fdebug+inode)
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_lehman_kernel(w, nw)
  !
  use constants, only : dp, eps12, cmplx_0
  use lattice,   only : ham
  use chi_internal, only: ek, ekq, occ_k, occ_kq, hk, hkq, chiff, chicc, chifc, chicf
  use IntRPA,    only : nFFidx, nCCidx, FFidx, CCidx
  !
  integer nw
  !
  complex(dp), dimension(nw) :: w
  !   w: input frequency
  !
  integer ibnd, jbnd
  integer ii, i1, i2
  real(dp) :: occ_diff
  complex(dp), dimension(nw) :: fact
  !
  complex(dp), dimension(nFFidx) :: UUf
  complex(dp), dimension(nCCidx) :: UUc
  !
  do ibnd=1, ham%norb
    do jbnd=1, ham%norb
      !
      occ_diff=occ_k(ibnd)-occ_kq(jbnd)
      !
      if (abs(occ_diff)>eps12) then
        !
        ! Substantial component
        fact(:)=occ_diff/(w(:)+ekq(jbnd)-ek(ibnd))
        !
        do ii=1, nFFidx
          !
          i1=FFidx(1, ii)
          i2=FFidx(2, ii)
          !
          UUf(ii)=conjg(hk(i1, ibnd))*hkq(i2, jbnd)
          !
        enddo
        !
        do ii=1, nCCidx
          !
          i1=CCidx(ii)
          !
          UUc(ii)=conjg(hk(i1, ibnd))*hkq(i1, jbnd)
          !
        enddo
        !
        ! call zgerc(m, n, alpha, x, incx, y, incy, a, lda)
        ! A=alpha * x * conjg(y') + A
        !
        do ii=1, nw
          !
          call zgerc(nFFidx, nFFidx, fact(ii), UUf, 1, UUf, 1, chiff(:, :, ii), nFFidx)
          !
          if (nCCidx>0) then
            !
            call zgerc(nFFidx, nCCidx, fact(ii), UUf, 1, UUc, 1, chifc(:, :, ii), nFFidx)
            call zgerc(nCCidx, nFFidx, fact(ii), UUc, 1, UUf, 1, chicf(:, :, ii), nCCidx)
            call zgerc(nCCidx, nCCidx, fact(ii), UUc, 1, UUc, 1, chicc(:, :, ii), nCCidx)
            !
          endif
          !
        enddo
        !
      endif
      !
    enddo ! jbnd
  enddo ! ibnd
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_lehman(chi0, w, nw, qv)
  !
  ! This subroutine calculates susceptibility matrix at qv
  !   of frequency w
  !
  use constants,    only : dp, cmplx_0, stdout
  use input,        only : beta
  use wanndata,     only : calc_hk
  use lattice,      only : nkirr, kvec, ham
  use chi_internal, only : hk, hkq, ek, ekq, occ_k, occ_kq, chiff, chicc, chifc, chicf
  use para,         only : para_merge_cmplx, distribute_calc, first_idx, last_idx
  use IntRPA,       only : nCCidx, nFFidx
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  real(dp),    dimension(3)  :: kq
  integer                    :: ik, iw
  !
  chiff=cmplx_0
  !
  if (nCCidx>0) then
    !
    chicc=cmplx_0
    chifc=cmplx_0
    chicf=cmplx_0
    !
  endif
  !
  call distribute_calc(nkirr)
  !
  do ik=first_idx, last_idx
    !
    !write(*, *) "ik: ", ik
    !write(*, *) "kvec: ", kvec(:, ik)
    call calc_hk(hk, ham, kvec(:, ik))
    !write(*, *) " hk: "
    !call show_matrix(hk, ham%norb, stdout)
    call calc_eig_occ_U(ek, occ_k, hk, beta)
    !
    kq=kvec(:, ik)+qv(:)
    call calc_hk(hkq, ham, kq)
    !write(*, *) " hkq: "
    !call show_matrix(hkq, ham%norb, stdout)
    call calc_eig_occ_U(ekq, occ_kq, hkq, beta)
    !
    call calc_chi_bare_matrix_lehman_kernel(w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chiff, nFFidx*nFFidx*nw)
  chiff=chiff/nkirr
  !
  if (nCCidx>0) then
    !
    call para_merge_cmplx(chicf, nCCidx*nFFidx*nw)
    call para_merge_cmplx(chifc, nFFidx*nCCidx*nw)
    call para_merge_cmplx(chicc, nCCidx*nCCidx*nw)
    !
    chifc=chifc/nkirr
    chicf=chicf/nkirr
    chicc=chicc/nkirr
    !
  endif
  !
  call calc_chi_trace_from_matrix(chi0, nw)
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_lehman_fast(chi0, w, nw, qv)
  !
  use constants,    only : dp, cmplx_0, stdout, twopi
  use wanndata,     only : calc_hk
  use lattice,      only : nkirr, kvec, ham
  use para,         only : para_merge_cmplx, first_idx, last_idx
  use IntRPA,       only : nCCidx, nFFidx
  use chi_internal, only : ek, ekq, Skq, occ_k, occ_kq, eig_full, occ_full, hk, hkq, Uk_local, Ukq_local, kq_idx, chiff, chifc, chicf, chicc
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  integer                    :: ik, ii, jj
  real(dp), dimension(3, last_idx-first_idx+1)  :: kq_local
  real(dp),    dimension(3)  :: DG
  real(dp)                   :: GdotTau
  complex(dp)                :: extra_phase
  complex(dp)                :: fact
  !
  chiff=cmplx_0
  !
  if (nCCidx>0) then
    !
    chicc=cmplx_0
    chifc=cmplx_0
    chicf=cmplx_0
    !
  endif
  !
  do ik=first_idx, last_idx
    !
    kq_local(:, ik-first_idx+1)=kvec(:, ik)+qv(:)
    !
  enddo
  ! Locate k+q elements in node 0 and populate
  call find_idx_in_kmesh(kq_local)
  call reorder_and_distribute_Ukq()
  !
  do ik=first_idx, last_idx
    !
    ek(:)=eig_full(:, ik)
    ekq(:)=eig_full(:, kq_idx(ik))
    !
    occ_k(:)=occ_full(:, ik)
    occ_kq(:)=occ_full(:, kq_idx(ik))
    !
    hk(:, :)=Uk_local(:, :, ik-first_idx+1)
    hkq(:, :)=Ukq_local(:, :, ik-first_idx+1)
    !
    ! Fix the phase due to site position
    DG(:)=kq_local(:, ik-first_idx+1)-kvec(:, kq_idx(ik))
    !
    do ii=1, ham%norb
      !
      GdotTau=SUM(DG(:)*ham%tau(:, ii))*twopi
      extra_phase=CMPLX(cos(GdotTau), sin(GdotTau), KIND=dp)
      hkq(ii, :)=conjg(extra_phase)*hkq(ii, :)
      !
    enddo
    !
    call calc_chi_bare_matrix_lehman_kernel(w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chiff, nFFidx*nFFidx*nw)
  chiff=chiff/nkirr
  !
  if (nCCidx>0) then
    !
    call para_merge_cmplx(chicf, nCCidx*nFFidx*nw)
    call para_merge_cmplx(chifc, nFFidx*nCCidx*nw)
    call para_merge_cmplx(chicc, nCCidx*nCCidx*nw)
    !
    chifc=chifc/nkirr
    chicf=chicf/nkirr
    chicc=chicc/nkirr
    !
  endif
  !
  call calc_chi_trace_from_matrix(chi0, nw)
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_GG_test(w, nw)
  !
  use constants,    only : dp, cmplx_1, cmplx_0, cmplx_i, stdout
  use lattice,      only : ham
  use input,        only : beta
  use chi_internal, only : hk, hkq, ek, ekq, chiff, chicc, chifc, chicf
  use IntRPA,       only : nFFidx, nCCidx, FFidx, CCidx
  use pade_sum,     only : npole, zp, eta
  !
  implicit none
  !
  integer                    :: nw
  complex(dp), dimension(nw) :: w
  !
  complex(dp), dimension(nw) :: fact
  complex(dp), dimension(ham%norb, ham%norb) :: gk, gkq
  complex(dp) :: z1, z2
  !
  integer iw, ii, i1, i2, jj, j1, j2, ipole
  !
  do ipole=1, npole
    !
    z1=zp(ipole)/beta*cmplx_i
    call calc_g0(gk, hk, z1, ham%norb, .false.)
    !
    do iw=1, nw
      !
      z2=z1+w(iw)
      call calc_g0(gkq, hkq, z2, ham%norb, .false.)
      !
      do ii=1, nFFidx
        !
        i1=FFidx(1, ii)
        i2=FFidx(2, ii)
        !
        do jj=1, nFFidx
          !
          j1=FFidx(1, jj)
          j2=FFidx(2, jj)
          !
          chiff(ii, jj, iw)=chiff(ii, jj, iw)+eta(ipole)*gk(i1, j1)*gkq(j2, i2)/beta
          !
        enddo ! jj
        !
      enddo   ! ii
      !
    enddo ! iw
    !
    z1=-zp(ipole)/beta*cmplx_i
    call calc_g0(gk, hk, z1, ham%norb, .false.)
    !
    do iw=1, nw
      !
      z2=z1+w(iw)
      call calc_g0(gkq, hkq, z2, ham%norb, .false.)
      !
      do ii=1, nFFidx
        !
        i1=FFidx(1, ii)
        i2=FFidx(2, ii)
        !
        do jj=1, nFFidx
          !
          j1=FFidx(1, jj)
          j2=FFidx(2, jj)
          !
          chiff(ii, jj, iw)=chiff(ii, jj, iw)+eta(ipole)*gk(i1, j1)*gkq(j2, i2)/beta
          !
        enddo ! jj
        !
      enddo   ! ii
      !
    enddo ! iw
    !
  enddo   ! ipole
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_GG_kernel(w, nw)
  !
  use constants,    only : dp, cmplx_1, cmplx_0, cmplx_i
  use lattice,      only : ham
  use input,        only : beta
  use chi_internal, only : hk, hkq, ek, ekq, chiff, chicc, chifc, chicf
  use IntRPA,       only : nFFidx, nCCidx, FFidx, CCidx
  use pade_sum,     only : npole, zp, eta
  !
  implicit none
  !
  integer                    :: nw
  complex(dp), dimension(nw) :: w
  !
  integer ibnd, jbnd
  integer ii, iw, i1, i2
  complex(dp), dimension(nw) :: fact
  complex(dp) :: z1, z2, g1, g2
  !
  complex(dp), dimension(nFFidx) :: UUf
  complex(dp), dimension(nCCidx) :: UUc
  !
  do ibnd=1, ham%norb
    do jbnd=1, ham%norb
      !
      fact(:)=cmplx_0
      !
      do ii=1, npole
        !
        z1=zp(ii)/beta*cmplx_i
        g1=cmplx_1/(z1-ek(ibnd))
        !
        do iw=1, nw
          z2=z1+w(iw)
          g2=cmplx_1/(z2-ekq(jbnd))
          fact(iw)=fact(iw)+eta(ii)*g1*g2
        enddo
        !
        z1=-zp(ii)/beta*cmplx_i
        g1=cmplx_1/(z1-ek(ibnd))
        !
        do iw=1, nw
          z2=z1+w(iw)
          g2=cmplx_1/(z2-ekq(jbnd))
          fact(iw)=fact(iw)+eta(ii)*g1*g2
        enddo
        !
      enddo
      !
      fact(:)=-fact(:)/beta
      !
      do ii=1, nFFidx
        !
        i1=FFidx(1, ii)
        i2=FFidx(2, ii)
        !
        UUf(ii)=conjg(hk(i1, ibnd))*hkq(i2, jbnd)
        !
      enddo
      !
      do ii=1, nCCidx
        !
        i1=CCidx(ii)
        !
        UUc(ii)=conjg(hk(i1, ibnd))*hkq(i1, ibnd)
        !
      enddo
      !write(*, *) "ibnd, jbnd:", ibnd, jbnd
      !write(*, '(2F14.9)') ek(ibnd), ekq(jbnd)
      !write(*, *) "UUf: "
      !write(*, '(18F9.4)') UUf
      !
      ! call zgerc(m, n, alpha, x, incx, y, incy, a, lda)
      ! A=alpha * x * conjg(y') + A
      !
      do ii=1, nw
        !
        call zgerc(nFFidx, nFFidx, fact(ii), UUf, 1, UUf, 1, chiff(:, :, ii), nFFidx)
        !
        !write(*, *) "chiFF : "
        !call show_matrix(chiff(:, :, 1), nFFidx, stdout)
        !
        if (nCCidx>0) then
          !
          call zgerc(nFFidx, nCCidx, fact(ii), UUf, 1, UUc, 1, chifc(:, :, ii), nFFidx)
          call zgerc(nCCidx, nFFidx, fact(ii), UUc, 1, UUf, 1, chicf(:, :, ii), nCCidx)
          call zgerc(nCCidx, nCCidx, fact(ii), UUc, 1, UUc, 1, chicc(:, :, ii), nCCidx)
          !
        endif
        !
      enddo
      !
    enddo ! jbnd
  enddo ! ibnd
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_GG(chi0, w, nw, qv)
  !
  ! This subroutine calculates susceptibility matrix at qv
  !   of frequency w
  !
  use constants,    only : dp, cmplx_0
  use input,        only : beta
  use wanndata,     only : calc_hk
  use lattice,      only : nkirr, kvec, ham
  use linalgwrap,   only : eigen
  use chi_internal, only : hk, hkq, ek, ekq, chiff, chicc, chifc, chicf
  use para,         only : para_merge_cmplx, distribute_calc, first_idx, last_idx
  use IntRPA,       only : nCCidx, nFFidx
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  real(dp),    dimension(3)  :: kq
  integer                    :: ik, iw
  !
  chiff=cmplx_0
  !
  if (nCCidx>0) then
    !
    chicc=cmplx_0
    chifc=cmplx_0
    chicf=cmplx_0
    !
  endif
  !
  call distribute_calc(nkirr)
  !
  do ik=first_idx, last_idx
    !
    !write(*, *) "inode, ik", inode, ik
    !write(*, *) "hk: "
    call calc_hk(hk, ham, kvec(:, ik))
    !call show_matrix(hk, ham%norb, stdout)
    call eigen(ek, hk, ham%norb)
    !write(*, *) "eig: "
    !write(*, '(10F9.4)') ek
    !write(*, *) "egv: "
    !call show_matrix(hk, ham%norb, stdout)
    !
    kq=kvec(:, ik)+qv(:)
    !write(*, *) "hkq:"
    call calc_hk(hkq, ham, kq)
    !call show_matrix(hkq, ham%norb, stdout)
    call eigen(ekq, hkq, ham%norb)
    !write(*, *) "eig: "
    !write(*, '(10F9.4)') ekq
    !write(*, *) "egv: "
    !call show_matrix(hkq, ham%norb, stdout)
    !
    call calc_chi_bare_matrix_GG_kernel(w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chiff, nFFidx*nFFidx*nw)
  chiff=chiff/nkirr
  !
  if (nCCidx>0) then
    !
    call para_merge_cmplx(chicf, nCCidx*nFFidx*nw)
    call para_merge_cmplx(chifc, nFFidx*nCCidx*nw)
    call para_merge_cmplx(chicc, nCCidx*nCCidx*nw)
    !
    chifc=chifc/nkirr
    chicf=chicf/nkirr
    chicc=chicc/nkirr
    !
  endif
  !
  call calc_chi_trace_from_matrix(chi0, nw)
  !
END SUBROUTINE

SUBROUTINE calc_chi_bare_matrix_GG_fast(chi0, w, nw, qv)
  !
  use constants,    only : dp, cmplx_0, twopi
  use wanndata,     only : calc_hk
  use lattice,      only : nkirr, kvec, ham
  use para,         only : para_merge_cmplx, first_idx, last_idx
  use IntRPA,       only : nCCidx, nFFidx
  use chi_internal, only : ek, ekq, eig_full, hk, hkq, Uk_local, Ukq_local, kq_idx, chiff, chifc, chicf, chicc
  !
  implicit none
  !
  integer nw
  !
  complex(dp), dimension(nw) :: chi0, w
  real(dp),    dimension(3)  :: qv
  !
  integer                    :: ik, ii, jj
  real(dp), dimension(3, last_idx-first_idx+1)  :: kq_local
  real(dp),    dimension(3)  :: DG
  real(dp)                   :: GdotTau
  complex(dp)                :: extra_phase
  !
  chiff=cmplx_0
  !
  if (nCCidx>0) then
    !
    chicc=cmplx_0
    chifc=cmplx_0
    chicf=cmplx_0
    !
  endif
  !
  do ik=first_idx, last_idx
    !
    kq_local(:, ik-first_idx+1)=kvec(:, ik)+qv(:)
    !
  enddo
  ! Locate k+q elements in node 0 and populate
  call find_idx_in_kmesh(kq_local)
  call reorder_and_distribute_Ukq()
  !
  do ik=first_idx, last_idx
    !
    ek(:)=eig_full(:, ik)
    ekq(:)=eig_full(:, kq_idx(ik))
    !
    hk(:, :)=Uk_local(:, :, ik-first_idx+1)
    hkq(:, :)=Ukq_local(:, :, ik-first_idx+1)
    !
    ! Fix the phase due to site position
    DG(:)=kq_local(:, ik-first_idx+1)-kvec(:, kq_idx(ik))
    !
    do ii=1, ham%norb
      !
      GdotTau=SUM(DG(:)*ham%tau(:, ii))*twopi
      extra_phase=CMPLX(cos(GdotTau), sin(GdotTau), KIND=dp)
      hkq(ii, :)=conjg(extra_phase)*hkq(ii, :)
      !
    enddo
    !
    call calc_chi_bare_matrix_GG_kernel(w, nw)
    !
  enddo
  !
  call para_merge_cmplx(chiff, nFFidx*nFFidx*nw)
  chiff=chiff/nkirr
  !
  if (nCCidx>0) then
    !
    call para_merge_cmplx(chicf, nCCidx*nFFidx*nw)
    call para_merge_cmplx(chifc, nFFidx*nCCidx*nw)
    call para_merge_cmplx(chicc, nCCidx*nCCidx*nw)
    !
    chifc=chifc/nkirr
    chicf=chicf/nkirr
    chicc=chicc/nkirr
    !
  endif
  !
  call calc_chi_trace_from_matrix(chi0, nw)
  !
END SUBROUTINE

SUBROUTINE calc_chi_trace_from_matrixFF_one(chi0)
  !
  use constants,    only : dp, cmplx_0, stdout
  use chi_internal, only : chiff
  use IntRPA,       only : nFFidx, FFidx
  !
  implicit none
  !
  complex(dp) :: chi0
  !
  integer ii, jj
  !
  chi0=cmplx_0
  !
  do ii=1, nFFidx
    if (FFidx(1, ii).ne.FFidx(2, ii)) cycle
    !
    do jj=1, nFFidx
      if (FFidx(1, jj).ne.FFidx(2, jj)) cycle
      !
      chi0=chi0+chiff(ii, jj, 1)
      !
    enddo
    !
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_chi_trace_from_matrixFF(chi0, nw)
  !
  use constants,    only : dp, cmplx_0, stdout
  use chi_internal, only : chiff
  use IntRPA,       only : nFFidx, FFidx
  !
  implicit none
  !
  integer nw
  complex(dp), dimension(nw) :: chi0
  !
  integer ii, jj
  !
  chi0=cmplx_0
  !
  do ii=1, nFFidx
    if (FFidx(1, ii).ne.FFidx(2, ii)) cycle
    !
    do jj=1, nFFidx
      if (FFidx(1, jj).ne.FFidx(2, jj)) cycle
      !
      chi0(:)=chi0(:)+chiff(ii, jj, :)
      !
    enddo
    !
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_chi_trace_from_matrix(chi0, nw)
  !
  use constants,    only : dp, cmplx_0, stdout
  use chi_internal, only : chiff, chicc, chifc, chicf
  use IntRPA,       only : nCCidx, nFFidx, FFidx, CCidx
  !
  implicit none
  !
  integer nw
  complex(dp), dimension(nw) :: chi0
  !
  integer ii, jj
  !
  chi0=cmplx_0
  !
  do ii=1, nFFidx
    if (FFidx(1, ii).ne.FFidx(2, ii)) cycle
    !
    do jj=1, nFFidx
      if (FFidx(1, jj).ne.FFidx(2, jj)) cycle
      !
      chi0(:)=chi0(:)+chiff(ii, jj, :)
      !
    enddo
    !
    do jj=1, nCCidx
      !
      chi0(:)=chi0(:)+chifc(ii, jj, :)
      chi0(:)=chi0(:)+chicf(jj, ii, :)
      !
    enddo
    !
  enddo
  !
  do ii=1, nCCidx
    do jj=1, nCCidx
      !
      chi0(:)=chi0(:)+chicc(ii, jj, :)
      !
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_eig_occ_U(eig, occ, hk, beta)
  !
  use constants,    only : dp, cmplx_1, cmplx_0
  use lattice,      only : ham
  use linalgwrap,   only : eigen
  !
  implicit none
  !
  real(dp), intent(in) :: beta
  complex(dp), dimension(ham%norb, ham%norb), intent(inout) :: hk
  !
  real(dp), dimension(ham%norb), intent(out) :: eig
  real(dp), dimension(ham%norb), intent(out) :: occ
  !
  real(dp) calc_occ
  !
  integer ii
  !
  call eigen(eig, hk, ham%norb)
  !
  do ii=1, ham%norb
    occ(ii)=calc_occ(eig(ii), beta)
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_Skq()
  !
  use constants,    only : dp, cmplx_1, cmplx_0
  use lattice,      only : ham
  use chi_internal, only : hk, hkq, Skq
  !
  implicit none
  !
  complex(dp), dimension(ham%norb, ham%norb) :: Stmp
  !
  integer ii, jj
  !
  ! calculates the structural function, 
  !  Skq(ibnd, jbnd)= | \sum_s conjg(egv_k(s, ibnd))*egv_kq(s, jbnd) |^2
  !   or in matrix form: Skq = | U_k^dagger * U_kq |^2
  !   !!! PLEASE BE AWARE THAT ibnd/jbnd order is relevant
  !       ibnd is the band index for k
  !    and jbnd is the band index for kq
  !
  !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call zgemm('C', 'N', ham%norb, ham%norb, ham%norb, cmplx_1, hk, ham%norb, hkq, ham%norb, cmplx_0, Stmp, ham%norb)
  !
  do ii=1, ham%norb
    do jj=1, ham%norb
      !Stmp(ii, jj)=sum(conjg(Uk(:, ii))*Ukq(:, jj))
      Skq(ii, jj)=abs(Stmp(ii, jj))**2
    enddo
  enddo
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

SUBROUTINE find_idx_in_kmesh(kq_local)
  !
  use constants,  only : dp, eps6
  use chi_internal, only : kq_idx
  use lattice,    only : nkirr, kvec, nk1, nk2, nk3
  use para,       only : first_idx, last_idx, para_merge_int, inode
  !
  implicit none
  !
  real(dp), dimension(3, last_idx-first_idx+1) :: kq_local
  !
  integer ii
  integer ik1, ik2, ik3
  !
  real(dp), dimension(3) :: kv
  !
  kq_idx(:)=0
  !
  do ii=1, last_idx-first_idx+1
    !
    ik1=mod(nint( (kq_local(1, ii)-floor(kq_local(1, ii)))*nk1 ), nk1)
    ik2=mod(nint( (kq_local(2, ii)-floor(kq_local(2, ii)))*nk2 ), nk2)
    ik3=mod(nint( (kq_local(3, ii)-floor(kq_local(3, ii)))*nk3 ), nk3)
    !
    kq_idx(ii+first_idx-1)=ik1*nk2*nk3+ik2*nk3+ik3+1
    !
    kv(1)=ik1*1.d0/nk1-kq_local(1, ii)
    kv(2)=ik2*1.d0/nk2-kq_local(2, ii)
    kv(3)=ik3*1.d0/nk3-kq_local(3, ii)
    !
    kv(:)=kv(:)-nint(kv(:))
    if (SUM(kv(:)*kv(:))>eps6) then
      write(*, *) " !!!! FATAL ERROR: K-MESH and QV inconsistent in FAST CALCULATION"
      stop
    endif
    !
  enddo
  !
  call para_merge_int(kq_idx, nkirr)
  !
END SUBROUTINE

SUBROUTINE reorder_and_distribute_Ukq()
  !
  use constants,   only : dp, cmplx_0
  use lattice,     only : ham, nkirr
  use chi_internal, only : Uk_full, Uk_local, Ukq_local, kq_idx
  use para,        only : first_idx, last_idx, para_distribute_cmplx, inode
  !
  implicit none
  !
  complex(dp), dimension(:, :, :), allocatable :: Ukq_full
  integer ii
  !
  if (inode.eq.0) then
    allocate(Ukq_full(ham%norb, ham%norb, nkirr))
    Ukq_full(:, :, :)=Uk_full(:, :, kq_idx(:))
  else
    allocate(Ukq_full(1, 1, 1))
  endif
  !
  call para_distribute_cmplx(Ukq_full, Ukq_local, ham%norb*ham%norb)
  !
  deallocate(Ukq_full)
  !
END SUBROUTINE

SUBROUTINE reorder_and_distribute_Gff()
  !
  use constants,   only : dp, cmplx_0
  use lattice,     only : ham, nkirr, ndimf
  use chi_internal, only : Gff_full, Gkqff_loc, kq_idx
  use para,        only : first_idx, last_idx, para_distribute_cmplx, inode
  use pade_sum,    only : npole
  !
  implicit none
  !
  complex(dp), dimension(:, :, :, :), allocatable :: Gkqff_full
  integer ii
  !
  if (inode.eq.0) then
    allocate(Gkqff_full(ndimf, ndimf, 2*npole, nkirr))
    Gkqff_full(:, :, :, :)=Gff_full(:, :, :, kq_idx(:))
  else
    allocate(Gkqff_full(1, 1, 1, 1))
  endif
  !
  call para_distribute_cmplx(Gkqff_full, Gkqff_loc, 2 * ndimf * ndimf * npole)
  !
  deallocate(Gkqff_full)
  !
END SUBROUTINE

SUBROUTINE reorder_and_distribute_Hff_Vfc_Ecc()
  !
  use constants,   only : dp, cmplx_0
  use lattice,     only : ham, nkirr, ndimf, ndimc
  use chi_internal, only : kq_idx, Hff_full, Vfc_full, Ecc_full, Hkqff_loc, Vkqfc_loc, Ecckq_loc
  use para,        only : first_idx, last_idx, para_distribute_cmplx, inode
  !
  implicit none
  !
  complex(dp), dimension(:, :, :), allocatable :: Mkq_full
  integer ii
  !
  ! Reorder Hff
  if (inode.eq.0) then
    allocate(Mkq_full(ndimf, ndimf, nkirr))
    Mkq_full(:, :, :)=Hff_full(:, :, kq_idx(:))
  else
    allocate(Mkq_full(1, 1, 1))
  endif
  !
  call para_distribute_cmplx(Mkq_full, Hkqff_loc, ndimf*ndimf)
  !
  deallocate(Mkq_full)
  !
  ! Reorder Vfc
  !
  if (inode.eq.0) then
    allocate(Mkq_full(ndimf, ndimc, nkirr))
    Mkq_full(:, :, :)=Vfc_full(:, :, kq_idx(:))
  else
    allocate(Mkq_full(1, 1, 1))
  endif
  !
  call para_distribute_cmplx(Mkq_full, Vkqfc_loc, ndimf*ndimc)
  !
  deallocate(Mkq_full)
  !
  Ecckq_loc(:, 1:last_idx-first_idx+1) = Ecc_full(:, kq_idx(first_idx:last_idx))
  !
END SUBROUTINE

SUBROUTINE show_matrix(mat, ndim, iou)
  !
  use constants,   only : dp
  !
  implicit none
  !
  integer :: ndim, iou
  complex(dp), dimension(ndim, ndim) :: mat
  !
  integer ii, jj
  !
  do ii=1, ndim
    !
    do jj=1, ndim
      !
      write(iou, '(2F9.4,1X)', advance='no') mat(ii, jj)
      !
    enddo
    !
    write(iou, *)
    !
  enddo
  !
END SUBROUTINE

SUBROUTINE show_rmatrix(mat, ndim, iou)
  !
  use constants,   only : dp
  !
  implicit none
  !
  integer :: ndim, iou
  real(dp), dimension(ndim, ndim) :: mat
  !
  integer ii, jj
  !
  do ii=1, ndim
    !
    do jj=1, ndim
      !
      write(iou, '(1F9.4,1X)', advance='no') mat(ii, jj)
      !
    enddo
    !
    write(iou, *)
    !
  enddo
  !
END SUBROUTINE