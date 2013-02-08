PROGRAM wannchi
  !
  use para
  use constants
  use wanndata, only: norb, finalize_wann
  use banddata, only: nbnd, nkx, nky, nkz, finalize_band
  !
  implicit none
  !
  integer nqpt
  integer iq
  real(dp) qvec(1:3)
  complex(dp), allocatable :: chi(:)
  character(len=80) seed
  !
  integer first_q, last_q
  !
  seed="wannier90"
  !
  CALL init_para
  !
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL read_input
  !
  CALL interpolate_bands
  !
  nqpt=nkx*nky
  allocate(chi(1:nqpt))
  !
  first_q=inode*nqpt/nnode+1
  last_q=(inode+1)*nqpt/nnode
  do iq=first_q, last_q
    !
    CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
    CALL compute_chi_diag(chi(iq), qvec)
    !
  enddo
  !
  CALL para_merge(chi, nqpt)
  !
  if (inode.eq.0) then
    do iq=1, nqpt
      CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
      !
      if (abs(qvec(1)).le.eps) write(*,*) ' '
      !
      write(*,'(3F12.8,2F22.12)') qvec, chi(iq)
    enddo
  endif
  !
  deallocate(chi)
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM
