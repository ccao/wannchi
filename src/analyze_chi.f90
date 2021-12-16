PROGRAM wannchi
  !
  use constants,only : stdout, dp, fout, cmplx_i, cmplx_0, cmplx_1, twopi
  use para,     only : init_para, inode, distribute_calc, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, mu, level, finalize_input, nnu, eps, nen, emesh
  use impurity, only : init_impurity, finalize_impurity, ismatsubara, nfreq, beta
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chi, chif, chic, w
  complex(dp) :: chi_diff
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
  allocate(chi(2*nfreq-nnu), chif(2*nfreq-nnu), chic(2*nfreq-nnu))
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
    !call analyze_chi_corr_trace1(chi, chif, chic, nnu, qvec(:, iq))
    call analyze_chi_corr_trace(chi, chif, nnu, qvec(:, iq))
    !
    if (inode.eq.0) then
      !
      write(fout, *) "# !!!!! ANALYZE OUTPUT "
      do ii=1, 2*nfreq-nnu
        chic(ii)=1.d0/chi(ii)-1.d0/chif(ii)
        write(fout, '(1F9.2,6G18.9)') ii-nfreq-0.5d0, chi(ii), chif(ii), chic(ii)
      enddo
      !
      !do ii=nfreq+1, 2*nfreq-nnu
      !  write(fout, '(1F19.12,2G18.9)') (ii-nfreq-0.5d0)*twopi/beta, chic(ii)+chic(2*nfreq+1-nnu-ii)
      !enddo
    endif
  enddo
  !
  deallocate(chif, chic)
  deallocate(chi)
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
