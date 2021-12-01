PROGRAM wannchi
  !
  use constants,only : cmplx_0, stdout, dp, fout
  use para,     only : init_para, inode, distribute_calc, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, mu, level, finalize_input, nnu
  use impurity, only : init_impurity, finalize_impurity, ismatsubara
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chi, chi0
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
  endif
  !
  allocate(chi(nnu), chi0(nnu))
  !
  if (inode.eq.0) then
    open(unit=fout, file='chi.dat')
    write(fout, '(A)') '#  susceptibility from wannchi'
  endif
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5)') 'Calculating qvec #', iq
    endif
    !
    call calc_chi_corr_trace(chi, chi0, nnu, qvec(:,iq))
    !
    if (inode.eq.0) then
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
  enddo
  !
  deallocate(chi, chi0)
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

