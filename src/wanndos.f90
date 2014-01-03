PROGRAM wanndos
  !
  use para
  use constants
  use wanndata, only: norb, finalize_wann
<<<<<<< HEAD
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
=======
  use banddata, only: egv, eig, nbnd, nkx, nky, nkz, finalize_band
  !
  implicit none
  !
  real(dp) emin, emax, de
  integer  nedos
  real(dp), allocatable :: dos(:), en(:)
  character(len=80) seed
  !
  integer first_e, last_e, idx_e, neig, ien
  !
>>>>>>> wanndos
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
<<<<<<< HEAD
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
=======
  deallocate(egv)
  !
  neig=nbnd*nkx*nky*nkz
  emin=min(eig)-0.5
  emax=max(eig)+0.5
  de=0.001
  sigma=0.01
  nedos=(emax-emin)/de
  allocate(dos(0:nedos))
  allocate(en(0:nedos))
  !
  do ien=0, nedos
    en(ien)=emin+de*ien
    dos(ien)=0.d0
  enddo
  !
  first_e=inode*neig/nnode+1
  last_e=(inode+1)*neig/nnode
  do idx_e=first_e, last_e
    !
    do ien=0, nedos
      if (abs(eig(idx_e)-en(ien))<0.02) then
        dos(ien)=dos(ien)+gaussian(eig(idx_e)-en(ien),sigma)
      endif
    enddo
    !
  enddo
  !
  CALL para_merge(chi, nqpt)
  !
  if (inode.eq.0) then
    do ien=0, nedos
      write(*,'(F12.3,2F22.12)') en(ien), dos(ien), SUM(dos(0:ien)*de)
    enddo
  endif
  !
  deallocate(dos)
>>>>>>> wanndos
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
<<<<<<< HEAD
=======
CONTAINS
  !
  function gaussian( x, sigma )
    !
    implicit none
    !
    real(dp) gaussian
    real(dp) x
    real(dp) sigma
    !
    gaussian=exp(-x*x/(sigma*sigma))/(sigma*1.77245385090552)
    !
  end function
  !
>>>>>>> wanndos
END PROGRAM
