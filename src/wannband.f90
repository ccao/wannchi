PROGRAM WannBand
  !
  use constants,    only : stdout, dp, fout
  use lattice,      only : ham, read_posfile, read_kmesh, read_impfile, read_sigfile, setup_mapping, fix_sigma_static, finalize_lattice_ham, finalize_lattice_kmesh, finalize_lattice_structure, print_impurity, interpolate_sigma_real, restore_lattice, nbath
  use wanndata,     only : read_ham, wannham_shift_ef, finalize_wann, calc_hk
  use para,         only : init_para, inode, distribute_calc, finalize_para, first_idx, last_idx, para_merge_real, para_sync_logical
  use input,        only : read_input, read_qpoints, mu, nqpt, qvec, nnu, nu, seed, finalize_input
  !
  implicit none
  !
  real(dp), dimension(:, :), allocatable :: Akw
  complex(dp), dimension(:, :), allocatable :: hk, hkw, gf
  complex(dp), dimension(:), allocatable :: sigma
  !
  logical isCorr
  !
  integer iq, ii
  !
  CALL init_para('WannBand')
  CALL read_input('wannband')
  !
  CALL read_ham(ham, seed)
  CALL wannham_shift_ef(ham, mu)
  !
  allocate(hk(ham%norb, ham%norb), gf(ham%norb, ham%norb))
  !
  if (inode.eq.0) then
    write(stdout, '(A,1F14.9,A)') "   # Fermi level shifted to ", mu, " eV"
    !
    inquire(file=trim(seed)//".impdef", exist=isCorr)
    !
  endif
  !
  call para_sync_logical(isCorr)
  !
  if (isCorr) then
    !
    CALL read_impfile(trim(seed)//".impdef")
    CALL setup_mapping
    !
    CALL read_sigfile(trim(seed)//".sig")
    CALL fix_sigma_static
    !
    allocate(sigma(nbath), hkw(ham%norb, ham%norb))
    !
    if (inode.eq.0) then
      write(stdout, *) " H0+Sinf: "
      call print_impurity(ham%hr(:, :, ham%r000))
    endif
    !
  endif
  !
  CALL read_posfile(trim(seed)//".pos")
  !
  !CALL read_kmesh("IBZKPT")
  CALL read_qpoints
  !
  allocate(Akw(nnu, nqpt))
  !
  call distribute_calc(nqpt)
  !
  do iq=first_idx, last_idx
    !
    if (inode.eq.0) then
      !
      write(stdout, '(A,1I5,A,1I5,A)', advance='no') " First node calculating ", iq, " out of ", (last_idx-first_idx+1), " q-points..."
      !
    endif
    !
    call calc_hk(hk, ham, qvec(:, iq))
    !
    do ii=1, nnu
      !
      if (isCorr) then
        !
        call interpolate_sigma_real(sigma, real(nu(ii)))
        call restore_lattice(hkw, sigma)
        !
        hkw(:, :)=hkw(:, :)+hk(:, :)
        !
        call calc_g0(gf, hkw, nu(ii), ham%norb, .false.)
        !
      else
        !
        call calc_g0(gf, hk, nu(ii), ham%norb, .false.)
        !
      endif
      !
      call calc_Akw(Akw(ii, iq), gf, ham%norb)
      !
    enddo ! ii
    !
    if (inode.eq.0) then
      !
      write(stdout, *) "done"
      !
    endif
    !
  enddo ! iq
  !
  if (inode.eq.0) then
    !
    write(stdout, *) " Generating output file..."
    !
  endif
  !
  call para_merge_real(Akw, nqpt*nnu)
  !
  if (inode.eq.0) then
    open(unit=fout, file="spectra.dat")
    !                    1...5...9.1...1...5...9.1...1...5...9.1...||1...5...9.1...||1...5...9.1...5...9.1.
    write(fout, '(A)') "# ================= qvec ================= || === w(i) === || ======== Akw ======== "
    !
    do iq=1, nqpt
      !
      do ii=1, nnu
      !
      write(fout, '(3F14.9,2X,1F14.9,2X,1G22.12)') qvec(:, iq), real(nu(ii)), Akw(ii, iq)
      !
      enddo
      !
    enddo
    !
    close(unit=fout)
    !
    write(stdout, *) " All done."
  endif
  !
  deallocate(hk, gf, Akw)
  !
  if (allocated(hkw)) deallocate(hkw)
  if (allocated(sigma)) deallocate(sigma)
  !
  CALL finalize_lattice_kmesh
  CALL finalize_lattice_structure
  CALL finalize_lattice_ham
  CALL finalize_input
  CALL finalize_para
  !
END PROGRAM
!
SUBROUTINE calc_Akw(Akw, gf, ndim)
  !
  use constants,  only : dp, twopi
  implicit none
  !
  integer ndim
  real(dp) :: Akw
  complex(dp), dimension(ndim, ndim) :: gf
  !
  integer ii
  !
  Akw=0.d0
  do ii=1, ndim
    Akw=Akw-aimag(gf(ii, ii))*2.d0/twopi
  enddo
  !
END SUBROUTINE