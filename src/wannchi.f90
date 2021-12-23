PROGRAM wannchi
  !
  use constants,only : stdout, dp, fout, cmplx_i, cmplx_0, eps4, twopi
  use para,     only : init_para, inode, distribute_calc, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, mu, level, finalize_input, nnu, eps, nen, emesh, beta, xq
  use impurity, only : init_impurity, finalize_impurity, ismatsubara, nfreq
  use chi_internal, only : init_chi_internal, finalize_chi_internal
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chiS, chiC, chi, chi0, w
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
  if (level>0) CALL init_impurity(beta, level)
  call init_chi_internal()
  if (level>1) CALL ham_fix_static()
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
  if (inode.eq.0) then
    open(unit=fout, file='chi0.dat')
    if (nnu>1) call output_header(fout)
    close(fout)
    if (level.eq.1) then
      open(unit=fout, file='chiS.dat')
      if (nnu>1) call output_header(fout)
      close(fout)
      open(unit=fout, file='chiC.dat')
      if (nnu>1) call output_header(fout)
      close(fout)
    endif
    if (level>1) then
      open(unit=fout, file='chiSigma.dat')
      if (nnu>1) call output_header(fout)
      close(fout)
    endif
  endif
  !
  allocate(chi0(nnu))
  if (level.eq.1) allocate(chiS(nnu), chiC(nnu))
  if (level>1) allocate(chi(nnu))
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5)') 'Calculating qvec #', iq
    endif
    !
    select case(level)
    case (0)
      call calc_chi_bare_trace(chi0, w, nnu, qvec(:, iq))
    case (1)
      call calc_chi_rpa_trace(chiS, chiC, chi0, w, nnu, qvec(:, iq))
    case (2)
      call calc_chi_bare_trace(chi0, w, nnu, qvec(:, iq))
    case (3)
      call calc_chi_corr_trace(chi, chi0, nnu, qvec(:,iq))
    end select
    !
    if (inode.eq.0) then
      open(unit=fout, file='chi0.dat', access='append')
      if (nnu>1) then
        call output_chi(chi0, fout, nqpt)
      else
        write(fout, '(3F14.9,2X,2G18.9)'), qvec(:, iq), chi0(1)
      endif
      close(fout)
      if (level.eq.1) then
        open(unit=fout, file='chiS.dat', access='append')
        if (nnu>1) then
          call output_chi(chiS, fout, nqpt)
        else
          write(fout, '(3F14.9,2X,2G18.9)'), qvec(:, iq), chiS(1)
        endif  
        close(fout)
        open(unit=fout, file='chiC.dat', access='append')
        if (nnu>1) then
          call output_chi(chiC, fout, nqpt)
        else
          write(fout, '(3F14.9,2X,2G18.9)'), qvec(:, iq), chiC(1)
        endif  
        close(fout)
      endif
      if (level>1) then
        open(unit=fout, file='chiSigma.dat', access='append')
        if (nnu>1) then
          call output_chi(chi, fout, nqpt)
        else
          write(fout, '(3F14.9,2X,2G18.9)'), qvec(:, iq), chi(1)
        endif  
        close(fout)
      endif
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
