PROGRAM wanndos
  !
  use para
  use constants
  use wanndata, only: norb, finalize_wann
  use banddata, only: eig, nbnd, nkx, nky, nkz, finalize_band, sigma
  use input, only: code
  !
  implicit none
  !
  real(dp) emin, emax, de, alpha
  integer  nedos
  real(dp) int_dos, dos, en
  real(dp), allocatable :: ene(:)
  character(len=80) seed
  !
  integer first_e, last_e, idx_e, neig, ien, ik
  !
  code=1
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
  neig=nbnd*nkx*nky*nkz
  !
  allocate(ene(1:neig))
  do ik=1, nkx*nky*nkz
    ene((ik-1)*nbnd+1:ik*nbnd)=eig(:,ik)
  enddo
  deallocate(eig)
  !
  emin=minval(ene)-0.5
  emax=maxval(ene)+0.5
  de=0.001
  nedos=(emax-emin)/de
  !
  do ien=0, nedos
    en=emin+de*ien
    dos=0.d0
    int_dos=0.d0
    first_e=inode*(neig/nnode)+1
    last_e=(inode+1)*(neig/nnode)
    do idx_e=first_e, last_e
      !
      alpha=(ene(idx_e)-en)/sigma
      dos=dos+2.d0*exp(-alpha*alpha)/sqrtpi
      int_dos=int_dos+erfc(alpha)
      !
    enddo
    !
    CALL para_merge(dos)
    CALL para_merge(int_dos)
    !
    dos=dos/(sigma*nkx*nky*nkz)
    int_dos=int_dos/(nkx*nky*nkz)
    !
    if (inode.eq.0) then
      write(*,'(F12.3,2F22.12)') en, dos, int_dos
    endif
  !
  enddo
  !
  deallocate(ene)
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM
