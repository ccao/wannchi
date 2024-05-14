PROGRAM wannchiRPA
  !
  use constants,    only : stdout, dp, fout
  use lattice,      only : read_posfile, ham
  use wanndata,     only : read_ham_dim, finalize_wann
  use para,         only : init_para, inode, distribute_calc, finalize_para
  use input,        only : read_input, read_qpoints, nqpt, qvec, ff_only, nnu, nu, seed, finalize_input
  use IntRPA,       only : read_RPA, finalize_RPA, calc_chiRPA
  use chi_internal, only : finalize_chi_internal, show_chi_diag, save_chi_matrix, read_chi_matrix_RPA, init_chi_matrix_RPA, chiff, chicc, chifc, chicf, chi0ff, chi0cc, chi0fc, chi0cf
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable  :: chi
  integer iq, ii
  !
  CALL init_para('WannChi RPA')
  CALL read_input('wannchi')
  !
  CALL read_ham_dim(ham, seed)
  !
  CALL read_posfile(trim(seed)//".pos")
  !
  CALL read_qpoints
  !
  CALL read_RPA
  !
  call init_chi_matrix_RPA(nnu)
  !
  allocate(chi(nnu))
  !
  open(unit=fout, file='chiRPAtr.dat')
  !
  !                    1...5...9.1...1...5...9.1...1...5...9.1...||1...5...9.1...1...5...9.1...||1...5...9.1...5...9.1.1...5...9.1...5...9.1.
  write(fout, '(A)') "# ================= qvec ================= || ========== w(i) ========== || =============== Tr[chiRPA] =============== "
  !
  do iq=1, nqpt
    !
    write(stdout, '(A,1I5,A,1I5,A)', advance='no') 'Calculating qvec #', iq, ' out of', nqpt, ' ...'
    !
    call read_chi_matrix_RPA((iq-1)*nnu, nnu, ff_only)
    !
    call calc_chiRPA(chiff, chicc, chifc, chicf, chi0ff, chi0cc, chi0fc, chi0cf, ff_only, nnu)
    call save_chi_matrix(iq*nnu, nnu, ff_only)
    !
    call calc_chi_trace_from_matrixFF(chi, nnu)
    !
    do ii=1, nnu
      !
      write(fout, '(3F14.9,2X,2F14.9,2X,2G22.12)') qvec(:, iq), nu(ii), chi(ii)
      !
    enddo
    !
    write(stdout, *) "done"
    !
  enddo
  !
  write(stdout, *) "All Done."
  !
  deallocate(chi)
  close(unit=fout)
  !
  CALL finalize_chi_internal
  CALL finalize_RPA
  CALL finalize_wann(ham, .true.)
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM

