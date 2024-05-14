PROGRAM PostChi
  !
  use constants,    only : stdout, dp, fout, cmplx_0
  use lattice,      only : read_posfile, ham, g2f_idx, partition, imp, read_impfile, setup_mapping
  use wanndata,     only : read_ham_dim, finalize_wann
  use para,         only : init_para, inode, distribute_calc, finalize_para
  use input,        only : read_input, read_qpoints, nqpt, qvec, ff_only, nnu, nu, seed, finalize_input
  use IntRPA,       only : read_RPA, finalize_RPA, nFFidx, FFidx
  use linalgwrap,   only : eigen
  use chi_internal, only : finalize_chi_internal, show_chi_diag, save_chi_matrix, read_chi_matrix, init_chi_matrix, chiff, chicc, chifc, chicf
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chi0
  real(dp), dimension(:),    allocatable :: chiEig
  integer, dimension(:),     allocatable :: spin
  integer iq, ii, jj, kk, i1, i2, j1, j2
  logical soc
  !
  CALL init_para('PostChi')
  CALL read_input('wannchi')
  !
  CALL read_ham_dim(ham, seed)
  CALL read_posfile(trim(seed)//".pos")
  CALL read_impfile(trim(seed)//".impdef")
  CALL setup_mapping
  !
  CALL read_qpoints
  !
  CALL read_RPA
  !
  ff_only=.true.
  !
  allocate(chiEig(nFFidx), chi0(nnu), spin(nFFidx))
  !
  call init_chi_matrix(nnu)
  !
  open(unit=fout, file='postchi.dat')
  !
  write(fout, '(A)') "# Up to 5 Principal eigenvalue and eigenstates. "
  !
  kk=partition(FFidx(1, 1))
  soc=(imp(kk)%ndim==(2*imp(kk)%lang+1))
  spin=0
  !
  if (soc) then
    !
    do ii=1, nFFidx
      i1=FFidx(1, ii)
      j1=partition(i1)
      kk=(imp(j1)%gidx(1)+imp(j1)%ndim/2-1)
      if (i1>kk) spin(ii)=spin(ii)+1
      i2=FFidx(2, ii)
      j2=partition(i2)
      kk=(imp(j2)%gidx(1)+imp(j2)%ndim/2-1)
      if (i2>kk) spin(ii)=spin(ii)+1
    enddo
    !
  endif
  !
  do iq=1, nqpt
    !
    write(fout, '(A,3F14.9)') "   #  Qvec: ", qvec(:, iq)
    !
    write(stdout, '(A,1I5,A,1I5,A)', advance='no') 'Calculating qvec #', iq, ' out of', nqpt, ' ...'
    !
    call read_chi_matrix((iq-1)*nnu, nnu, ff_only)
    !
    if (soc) then
      do ii=1, nFFidx
        do jj=1, nFFidx
          !
          if (MOD(spin(ii)+spin(jj), 2)==1) chiff(ii, jj, :)=cmplx_0
          !
        enddo
      enddo
    endif
    !
    call calc_chi_trace_from_matrixFF(chi0, nnu)
    !
    do ii=1, nnu
      !
      write(fout, '(A,2F14.9,A,2G22.12)') "   #    w= ", nu(ii), " Tr[chi0] : ", chi0(ii)
      !
      call eigen(chiEig, chiff(:, :, ii), nFFidx)
      !
      if (nFFidx>4) then
        do jj=0, 4
          write(fout, '(A,1I4,A,1G22.12)') "            chiEig ", jj+1, " : ", chiEig(nFFidx-jj)
          do kk=1, nFFidx
            if (abs(chiff(kk, nFFidx-jj, ii))>0.01) write(fout, '(15X,2I5,2F14.9)') FFidx(1, kk), FFidx(2, kk), chiff(kk, nFFidx-jj, ii)
          enddo
        enddo
      else
        do jj=0, nFFidx-1
          write(fout, '(A,1I4,A,1G22.12)') "            chiEig ", jj+1, " : ", chiEig(nFFidx-jj)
          do kk=1, nFFidx
            write(fout, '(15X,2I5,2F14.9)') FFidx(1, kk), FFidx(2, kk), chiff(kk, nFFidx-jj, ii)
          enddo
        enddo
      endif
      !
      write(fout, '(A)') "   #          "
      !
    enddo
    !
    write(fout, '(A)') "   #          "
    !
    write(stdout, *) "done"
    !
  enddo
  !
  write(stdout, *) "All Done."
  !
  close(unit=fout)
  !
  deallocate(chiEig)
  !
  CALL finalize_chi_internal
  CALL finalize_RPA
  CALL finalize_wann(ham, .true.)
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM

