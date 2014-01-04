PROGRAM wannchi
  !
  use para
  use constants
  use wanndata, only: norb, finalize_wann
  use banddata, only: nbnd, nkx, nky, nkz, finalize_band
  use input, only: code, mode, nqseg, nqbnd, get_qvec, lrpa
  use chidata, only: finalize_chi
  !
  implicit none
  !
  integer nqpt
  integer iq
  real(dp) qvec(1:3)
  complex(dp), allocatable :: chi(:)
  character(len=80) seed
  !
  integer first_q, last_q, dnq
  !
  code=0
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
  if (lrpa) then
    CALL compute_U_mat
  endif
  !
  CALL interpolate_bands
  !
  if (mode.eq.0) then
    nqpt=nkx*nky
  else ! (mode.eq.1)
    nqpt=nqbnd*nqseg+1
  endif

  first_q=inode*nqpt/nnode+1
  last_q=(inode+1)*nqpt/nnode

  allocate(chi(1:nqpt))
  !
  dnq=(last_q-first_q)/100
  if (dnq<1) dnq=2
  !
  do iq=first_q, last_q
    !
    if(mode.eq.0) then
      CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
    else ! (mode.eq.1)
      CALL get_qvec(qvec, iq)
    endif
    CALL compute_chi_diag(chi(iq), qvec)
    !
    if (mod(iq-first_q,dnq).eq.0) then
      if (inode.eq.0) then
        write(*,*) " #... Percentage done: ", (iq-first_q)*100/(last_q-first_q), "%"
      endif
    endif
    !
  enddo
  !
  CALL para_merge(chi, nqpt)
  !
  if (inode.eq.0) then
    do iq=1, nqpt
      if(mode.eq.0) then
        CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
      else
        CALL get_qvec(qvec, iq)
      endif
      !
      if (abs(qvec(1)).le.eps) write(*,*) ' '
      !
      write(*,'(3F12.8,2F22.12)') qvec, chi(iq)
    enddo
  endif
  !
  deallocate(chi)
  !
  CALL finalize_chi
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM
