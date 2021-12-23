MODULE input
  !
  use constants
  !
  implicit none
  !
  real(dp), dimension(3, 3) :: acell
  ! Lattice vectors acell(3, 3)
  real(dp), dimension(3, 3) :: bcell
  ! Reciprocal lattice vectors bcell(3,3)
  integer nqpt
  ! Number of qpts
  real(dp), dimension(:, :), allocatable :: qvec
  ! qvectors: qvec(3, iq)
  real(dp) mu, beta, eps
  !
  integer nnu
  ! Number of matsubara frequencies to be calculated
  integer nexact
  ! Number of exact Matsubara frequency GFs to be calculated
  !
  integer  nkx, nky, nkz
  ! Kmesh in interpolation
  integer  nen
  real(dp) :: emin, emax
  real(dp), dimension(:), allocatable :: emesh
  !
  !These are used to output the header
  integer nqsec
  real(dp), dimension(:), allocatable :: xq
  character(len=5), dimension(:), allocatable :: xlabel
  !
  integer level
  ! Calculation level
  integer mode
  ! Calculation mode
  !
  character(len=80) seed
  !
  namelist /CONTROL/ level, mode, nkx, nky, nkz, nnu, nexact, nen, emin, emax, eps
  namelist /SYSTEM/ seed, beta, mu
  !
CONTAINS
  !
 SUBROUTINE read_input(codename)
  !
  !
  use constants, only : dp, eps4, fin, fout
  use para,      only : inode, para_sync_int, para_sync_real
  use linalgwrap, only : invmat
  !
  implicit none
  !
  character(*) codename
  !
  character(len=80) line
  integer nq_per_sec, iqsec, iq
  integer iqx, iqy, iqz
  real(dp), allocatable :: qbnd_vec(:, :)
  real(dp), dimension(5) :: tt_real
  integer, dimension(8) :: tt_int
  !character dir
  ! INPUT FILE:
  ! &SYSTEM
  !   seed='wannier90',
  !   beta=2000.d0,
  !   mu=0.d0,
  ! /
  ! &CONTROL
  !   level=2
  !   mode=1
  !   nkx=48
  !   nky=48
  !   nkz=48
  !   nnu=64
  !   nexact=4000
  !   nen=5001
  !   emin=-3.d0
  !   emax=2.d0
  !   eps=0.001
  ! /
  ! CELL
  !   1.d0 0.d0 0.d0
  !   0.d0 1.d0 0.d0
  !   0.d0 0.d0 1.d0
  ! QPOINTS
  !   0.0 0.0 0.0
  !
  ! DEFAULT:
  level=0
  ! Noninteracting caculation
  ! 1: RPA
  ! 2: s_oo (static fix)
  ! 3: dynamic (self energy correction)
  mode=0
  ! Single point calculation
  nkx=1
  nky=1
  nkz=1
  ! Gamma Point
  nnu=1
  ! Calculate only single frequency
  nexact=0
  !
  nen=1
  ! Single frequency
  emin=0.d0
  emax=0.d0
  ! Default at Fermi level
  seed='wannier90'
  ! Default seed name
  beta=1e7
  ! Very low T calculation (~0K)
  mu=0.d0
  ! Default Fermi level 
  eps=eps4
  ! Default small imaginary part
  if (inode.eq.0) then
    open(unit=fin, file=trim(codename)//".inp")
    !
    read(nml=SYSTEM, unit=fin)
    read(nml=CONTROL, unit=fin)
    !
    tt_int(:)=(/level, mode, nkx, nky, nkz, nnu, nexact, nen/);
    tt_real(:)=(/emin, emax, eps, beta, mu/);
    !
    read(fin, *) line
    !
    if (trim(line)/='CELL') then
      write(*, *) "!!! FATAL ERROR: Input format incorrect, CELL?"
      stop
    endif
    do iq=1, 3
      read(fin, *) acell(:, iq)
    enddo
    read(fin, *) line
    if (trim(line)/='QPOINTS') then
      write(*, *) "!!! FATAL ERROR: Input format incorrect, QPOINT?"
      stop
    endif
    !
    if (level==0) write(stdout, '(A)', advance='no') "    Uncorrelated calculation "
    if (level==1) write(stdout, '(A)', advance='no') "    RPA calculation "
    if (level==2) write(stdout, '(A)', advance='no') "    Correlated static limit calculation "
    if (level==3) write(stdout, '(A)', advance='no') "    Correlated dynamic calculation "
    !
  endif
  !
  call para_sync_int(tt_int, 8)
  level=tt_int(1)
  mode=tt_int(2)
  nkx=tt_int(3)
  nky=tt_int(4)
  nkz=tt_int(5)
  nnu=tt_int(6)
  nexact=tt_int(7)
  nen=tt_int(8)
  !
  call para_sync_real(tt_real, 5)
  emin=tt_real(1)
  emax=tt_real(2)
  eps=tt_real(3)
  beta=tt_real(4)
  mu=tt_real(5)
  !
  allocate(emesh(nen))
  !
  if (nen>1) then
    do iq=1, nen
      emesh(iq)=emin+(emax-emin)*(iq-1)/(nen-1)
    enddo
  else
    emesh(1)=emin
  endif
  !
  call para_sync_real(acell, 9)
  bcell=acell
  call invmat(bcell, 3)
  !
  if (beta<0) beta=1.d7
  if (eps>eps4.or.eps<eps9) eps=eps4
  !
  select case (mode)
    case (0)
      nqpt=1
      allocate(qvec(1:3, 1:1))
      if (inode.eq.0) then
        write(stdout, '(A)') "of a single point"
        read(fin, *) qvec(:, 1)
      endif
      call para_sync_real(qvec, 3*nqpt)
    case (1)
      if (inode.eq.0) then
        write(stdout, '(A)') "of band structure"
        read(fin, *) tt_int(1:2)
      endif
      call para_sync_int(tt_int, 2)
      nqsec=tt_int(1)
      nq_per_sec=tt_int(2)
      nqpt=(nqsec-1)*nq_per_sec+1
      allocate(qbnd_vec(1:3, 1:nqsec))
      allocate(qvec(1:3, nqpt))
      allocate(xq(nqsec))
      allocate(xlabel(nqsec))
      !
      if (inode.eq.0) then
        do iqsec=1, nqsec
          read(fin, *) qbnd_vec(:, iqsec), xlabel(iqsec)
          if (iqsec==1) then
            xq(iqsec)=0.d0
          else
            tt_real(1:3)=qbnd_vec(:, iqsec)-qbnd_vec(:,iqsec-1)
            do iqx=1, 3
              tt_real(3+iqx)=sum(bcell(iqx,:)*tt_real(1:3))
            enddo
            xq(iqsec)=xq(iqsec-1)+sqrt(sum(tt_real(4:6)**2))
          endif
        enddo
        !
        do iqsec=1, nqsec-1
          do iq=1, nq_per_sec
            qvec(:, (iqsec-1)*nq_per_sec+iq)=qbnd_vec(:, iqsec)+(qbnd_vec(:, iqsec+1)-qbnd_vec(:, iqsec))*(iq-1)/nq_per_sec
          enddo
        enddo
        qvec(:, nqpt)=qbnd_vec(:, nqsec)
        deallocate(qbnd_vec)
      endif
      !
      call para_sync_real(qvec, 3*nqpt)
      call para_sync_real(xq, nqsec)
    case (2)
      if (inode.eq.0) then
        write(stdout, '(A)') "of a plane"
      endif
      nqpt=nkx*nky
      allocate(qvec(1:3, 1:nqpt))
      do iqx=1, nkx
        do iqy=1, nky
          iq=(iqy-1)*nkx+iqx
          qvec(1, iq)=(iqx-1.d0)/nkx
          qvec(2, iq)=(iqy-1.d0)/nky
          qvec(3, iq)=0.d0
        enddo
      enddo
    case (3)
      if (inode.eq.0) then
        write(stdout, '(A)') "of bulk"
      endif
      nqpt=nkx*nky*nkz
      allocate(qvec(1:3, 1:nqpt))
      do iqx=1, nkx
        do iqy=1, nky
          do iqz=1, nkz
            iq=(iqz-1)*nkx*nky+(iqy-1)*nkx+iqx
            qvec(1, iq)=(iqx-1.d0)/nkx
            qvec(2, iq)=(iqy-1.d0)/nky
            qvec(3, iq)=(iqz-1.d0)/nkz
          enddo
        enddo
      enddo
  end select
  !
  if (inode.eq.0) close(unit=fin)
  !
 END SUBROUTINE
  !
 SUBROUTINE finalize_input
  !
  implicit none
  !
  if (allocated(qvec)) deallocate(qvec)
  if (allocated(emesh)) deallocate(emesh)
  if (allocated(xq)) deallocate(xq)
  if (allocated(xlabel)) deallocate(xlabel)
  !
 END SUBROUTINE
  !
END MODULE
