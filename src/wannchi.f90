PROGRAM wannchi
  !
  use constants,    only : stdout, dp, fout
  use lattice,      only : read_posfile, read_kmesh, ham
  use wanndata,     only : read_ham, wannham_shift_ef, finalize_wann
  use para,         only : init_para, inode, distribute_calc, finalize_para
  use input,        only : read_input, read_qpoints, mu, nqpt, qvec, use_lehman, fast_calc, trace_only, ff_only, nnu, nu, seed, finalize_input
  use IntRPA,       only : read_RPA
  use pade_sum,     only : print_pade
  !
  implicit none
  !
  complex(dp), dimension(:), allocatable :: chi0
  integer iq, ii
  !
  CALL init_para('WannChi')
  CALL read_input('wannchi')
  !
  CALL read_ham(ham, seed)
  CALL wannham_shift_ef(ham, mu)
  !
  CALL read_kmesh("IBZKPT")
  !
  CALL read_qpoints
  !
  if (.not.trace_only) CALL read_RPA
  !
  if (inode.eq.0) then
    write(stdout, '(A,1F14.9,A)') "   # Fermi level shifted to ", mu, " eV"
  endif
  !
  if (use_lehman) then
    !
    call prepare_lehman
    !
  else
    !
    call prepare_GG
    !
    call print_pade
    !
  endif
  !
  allocate(chi0(nnu))
  !
  if (inode.eq.0) then
    !
    open(unit=fout, file='chi0tr.dat')
    !
    !                    1...5...9.1...1...5...9.1...1...5...9.1...||1...5...9.1...1...5...9.1...||1...5...9.1...5...9.1.1...5...9.1...5...9.1.
    write(fout, '(A)') "# ================= qvec ================= || ========== w(i) ========== || ================ Tr[chi0] ================ "
    !
  endif
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5)') 'Calculating qvec #', iq
    endif
    !
    if (use_lehman) then
      !
      if (trace_only) then
        !
        if (fast_calc) then
          call calc_chi_bare_trace_lehman_fast(chi0, nu, nnu, qvec(:, iq))
        else
          call calc_chi_bare_trace_lehman(chi0, nu, nnu, qvec(:, iq))
        endif
        !
      else
        ! Full matrix
        if (fast_calc) then
          call calc_chi_bare_matrix_lehman_fast(chi0, nu, nnu, qvec(:, iq))
        else
          call calc_chi_bare_matrix_lehman(chi0, nu, nnu, qvec(:, iq))
        endif
        !
      endif
      !
    else
      ! Use G*G
      if (fast_calc) then
        call calc_chi_bare_matrix_GG_fast(chi0, nu, nnu, qvec(:, iq))
      else
        call calc_chi_bare_matrix_GG(chi0, nu, nnu, qvec(:, iq))
      endif
      !
    endif
    !
    if (inode.eq.0) then
      !
      do ii=1, nnu
        write(fout, '(3F14.9,2X,2F14.9,2X,2G22.12)') qvec(:, iq), nu(ii), chi0(ii)
      enddo
      !
      if (.not.trace_only) then
        ! Full matrix to binary files
        !
        if (ff_only) then
          !
          !
        else
          !
          !
        endif
        !
      endif
      !
    endif
  enddo
  !
  deallocate(chi0)
  !
  if (inode.eq.0) close(unit=fout)
  !
  CALL finalize_wann(ham, .true.)
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM
