PROGRAM wannchi
  !
  use constants,only : cmplx_0, stdout, dp, fout
  use para,     only : init_para, inode, distribute_k, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann
  use banddata, only : nbnd, init_band, finalize_band
  use input,    only : read_input, seed, qvec, nqpt, nkpt
  use chidata,  only : finalize_chi
  !
  implicit none
  !
  integer iq
  integer dnq
  complex(dp), allocatable :: chi(:)
  !
  CALL init_para
  CALL read_input
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL init_band
  !
  CALL distribute_k(nkpt)
  !
  CALL interpolate_bands
  !
  allocate(chi(1:nqpt))
  chi=cmplx_0
  !
  open(unit=fout, file="wannchi.dat")
  !
  dnq=nqpt/100
  if (dnq<2) dnq=2
  !
  do iq=1, nqpt
    !
    CALL compute_chi_bare_diag(chi(iq), qvec(:, iq))
    !
    if ((nqpt.ne.1).and.(mod(iq-1, dnq).eq.0)) then
      if (inode.eq.0) then
        write(stdout, *) " #... Percentage done: ", (iq-1)*100/nqpt, "%"
      endif
    endif
    !
    CALL output_chi(qvec(:, iq), chi(iq))
    !
  enddo
  !
  deallocate(chi)
  !
  close(unit=fout)
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
