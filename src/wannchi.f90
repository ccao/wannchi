PROGRAM wannchi
  !
  use constants,only : stdout, dp, fout, cmplx_i, cmplx_0
  use para,     only : init_para, inode, distribute_calc, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, mu, level, finalize_input, nnu, eps, nen, emesh
  use impurity, only : init_impurity, finalize_impurity, ismatsubara, nfreq
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
  CALL ham_shift_ef(mu)
  !
  if (level>0) then
    CALL init_impurity()
    CALL ham_fix_static()
  endif
  !
  if (ismatsubara) then
    call check_beta()
  else
    nen=nnu
    allocate(emesh(nnu), w(nnu))
    do ii=1, nnu
      emesh(ii)=(ii-1.d0)*eps
      if (ii>1) then
        w(ii)=emesh(ii)+eps*0.1*cmplx_i
      else
        w(ii)=cmplx_0
      endif
    enddo
  endif
  !
  allocate(chi0(nnu))
  if (level>1) allocate(chi(nnu))
  !allocate(chi(2*nfreq-nnu), chi0(2*nfreq-nnu))
  !
  if (inode.eq.0) then
    open(unit=fout, file='chi.dat')
    if (level<2) call output_header(fout)
  endif
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5)') 'Calculating qvec #', iq
    endif
    !
    if (level>1) then
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
      !write(fout, *) "# !!!!! ANALYZE OUTPUT "
      !do ii=1, 2*nfreq-nnu
      !  chi_diff=cmplx_1/chi(ii)-cmplx_1/chi0(ii)
      !  write(fout, '(1F9.2,6G18.9)') ii-nfreq-0.5d0, chi(ii), chi0(ii), chi_diff
      !enddo
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

SUBROUTINE check_beta()
  !
  use constants, only : eps4, stdout
  use impurity,  only : beta
  use input,     only : beta_in
  !
  implicit none
  !
  if (abs(beta-beta_in)>eps4) then
    write(stdout, '(A)') '!!! WARNNING: input beta is not consistent with self-energy!'
    write(stdout, '(A,2F9.4)') '        The input and sig file beta are:', beta_in, beta
    write(stdout, '(A)') '      We shall continue using the input self-energy'
    write(stdout, '(A)') '      But you must know what you are doing!'
  endif
  !
  beta=beta_in
  !
END SUBROUTINE

