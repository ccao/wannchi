PROGRAM wannchi
  !
  use constants
<<<<<<< HEAD
  use wanndata, only: norb, finalize_wann
  use banddata, only: nbnd, nkx, nky, nkz, finalize_band
  use input, only: code, mode, nqseg, nqbnd, get_qvec, lrpa
  use chidata, only: finalize_chi
=======
  use para,     only : init_para
  use wanndata, only : read_ham, norb
  use banddata, only : nkx, nky, nkz, nbnd
  use input,    only : read_input, seed, nqx, nqy, nqz
>>>>>>> New modulized version
  !
  implicit none
  !
  integer iq
  integer first_q, last_q, dnq
  complex(dp), allocatable :: chi(:)
<<<<<<< HEAD
  character(len=80) seed
  !
  integer first_q, last_q, dnq
  !
  code=0
  seed="wannier90"
=======
>>>>>>> New modulized version
  !
  CALL init_para
  CALL read_input
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL init_band
  !
  if (lrpa) then
    CALL compute_U_mat
  endif
  !
<<<<<<< HEAD
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
=======
  first_q=inode*nqpt/nnode+1
  last_q=(inode+1)*nqpt/nnode
  !
  allocate(chi(1:nqpt))
  !
  dnq=(last_q-first_q)/100
  if (dnq<2) dnq=2
  !
  do iq=first_q, last_q
    !
    CALL compute_chi_bare_diag(chi(iq), qvec(:, iq))
    !
    if (mod(iq-first_q, dnq).eq.0) then
      if (inode.eq.0) then
        write(stdout, *) " #... Percentage done: ", (iq-first_q)*100/(last_q-fist_q), "%"
      endif
    endif
>>>>>>> New modulized version
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
<<<<<<< HEAD
      if(mode.eq.0) then
        CALL map_to_kpt(qvec, iq, nkx, nky, nkz)
      else
        CALL get_qvec(qvec, iq)
      endif
      !
      if (abs(qvec(1)).le.eps) write(*,*) ' '
      !
      write(*,'(3F12.8,2F22.12)') qvec, chi(iq)
=======
      write(stdout, '(3F12.8,2F22.12)') qvec(:, iq), chi(iq)
>>>>>>> New modulized version
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
