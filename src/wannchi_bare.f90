PROGRAM wannchi
  !
  use constants,only : stdout, dp, fout, cmplx_i, cmplx_0, eps4, twopi
  use para,     only : init_para, inode, distribute_calc, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, mu, level, finalize_input, nnu, eps, nen, emesh, beta
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chi, chi0, w
  integer iq, ii
  !
  CALL init_para('wannchi')
  CALL read_input('wannchi')
  CALL read_ham(seed)
  !
  if (inode.eq.0) then
    write(stdout, '(A,1F14.9,A)') "   # Fermi level shifted to ", mu, " eV"
  endif
  !
  CALL ham_shift_ef(mu)
  !
  if (level>0) then
    CALL init_impurity(beta)
    if (level>1) then
      CALL ham_fix_static()
    else
      CALL 
  endif
  !
  if (inode.eq.0) then
    open(unit=fout, file='chi.dat')
    if (level<3) call output_header(fout)
  endif
  !
  if (ismatsubara) then
    allocate(w(nnu))
    do ii=1, nnu
      w(ii)=(ii-1.d0)*twopi/beta*cmplx_i
    enddo
  else
    nnu=nen
    allocate(w(nnu))
    w(:)=emesh(:)+eps*cmplx_i
  endif
  !
  allocate(chi0(nnu))
  if (level>1) allocate(chi(nnu))
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5)') 'Calculating qvec #', iq
    endif
    !
    if (level>2) then
      call calc_chi_corr_trace(chi, chi0, nnu, qvec(:,iq))
      !call analyze_chi_corr_trace(chi, chi0, nnu, qvec(:, iq))
    else
      call calc_chi_bare_trace(chi0, w, nnu, qvec(:, iq))
    endif
    !
    if (inode.eq.0) then
      if (level<2) then
        call output_chi(chi0, fout, nqpt)
      else
        write(fout, '(A,3F9.4)') '#   q: ', qvec(:, iq)
        write(fout, '(A)') '  # Non-interacting chi0:'
        do ii=1, nnu
          write(fout, '(2F14.9,2X)', advance='no') chi0(ii)
        enddo
        write(fout, *)
        write(fout, '(A)') '  # Interacting chi: '
        do ii=1, nnu
          write(fout, '(2F14.9,2X)', advance='no') chi(ii)
        enddo
        write(fout, *)
      endif
      !
    endif
  enddo
  !
  deallocate(chi0)
  if (level>1) deallocate(chi)
  !
  if (inode.eq.0) close(unit=fout)
  !
  CALL finalize_wann
  CALL finalize_impurity
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM
