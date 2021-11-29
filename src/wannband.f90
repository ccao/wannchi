PROGRAM wannband
  !
  use constants,only : cmplx_0, stdout, dp, fout, cmplx_i, twopi
  use para,     only : init_para, inode, distribute_calc, finalize_para, first_idx, last_idx, para_merge_real
  use wanndata, only : read_ham, norb, finalize_wann, ham_shift_ef
  use input,    only : read_input, seed, qvec, nqpt, emesh, nen, level, eps, finalize_input
  use impurity, only : sinf, fulldim, restore_lattice, basis_map, init_impurity
  !
  implicit none
  !
  real(dp), dimension(:), allocatable :: spec
  complex(dp), dimension(:, :), allocatable :: hk
  complex(dp), dimension(:, :), allocatable :: gf
  complex(dp) :: w
  !
  integer iq, ii, jj
  integer dnq
  !
  CALL init_para('wannband')
  CALL read_input('wannband')
  CALL read_ham(seed)
  !
  if (level==2) then
    CALL init_impurity()
    CALL ham_fix_static()
  endif
  !
  CALL distribute_calc(nen)
  !
  allocate(hk(norb, norb))
  allocate(gf(norb, norb))
  allocate(spec(nen))
  !
  if (inode.eq.0) then
    open(unit=fout, file='spectra.dat')
    call output_header(fout)
  endif
  !
  dnq=nqpt/100
  if (dnq<2) dnq=2
  !
  do iq=1, nqpt
    !
    if (inode.eq.0) then
      if (mod(iq-1, dnq).eq.0) then
        write(stdout, *) "  ... Percentage done: ", (iq-1)*100/nqpt, "%"
      endif
    endif
    !
    spec(:)=0.d0
    call calc_hk(hk, qvec(:, iq))
    !
    do ii=first_idx, last_idx
      !
      w=emesh(ii)+eps*cmplx_i
      !
      if (level==2) then
         call calc_corr_realgf(gf, hk, w, .false.)
      else
        call calc_g0(gf, hk, w, norb, .false.)
      endif
      !
      do jj=1, norb
        spec(ii)=spec(ii)-aimag(gf(jj, jj))/twopi
      enddo
      !
    enddo ! ii
    !
    call para_merge_real(spec, nen)
    !
    if (inode.eq.0) then
      call output_spectral(spec, emesh, nen, fout, nqpt)
    endif
    !
  enddo ! iq
  !
  if (inode.eq.0) close(unit=fout)
  !
  deallocate(hk, gf, spec)
  !
  CALL finalize_wann
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM