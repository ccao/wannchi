MODULE input
  !
  use constants
  !
  implicit none
  !
  integer nqpt
  real(dp), dimension(3, 3) :: acell
  ! Lattice vectors acell(3, 3)
  real(dp), dimension(3, 3) :: bcell
  ! Reciprocal lattice vectors bcell(3,3)
  real(dp), dimension(:, :), allocatable :: qvec
  ! qvectors: qvec(3, iq)
  real(dp) mu, beta, omega, eps
  integer  nkx, nky, nkz
  integer  nen
  real(dp), dimension(:), allocatable :: emesh
  !
  !These are used to output the header
  integer nqsec
  real(dp), dimension(:), allocatable :: xq
  character(len=5), dimension(:), allocatable :: xlabel
  !
  integer level
  integer mode
  !
  character(len=80) seed
  !
CONTAINS
  !
 SUBROUTINE read_input(codename)
  !
  !************ INPUT FILE *************
  !** file name: wannchi.inp
  !line 1: seed name       ! seed name
  !line 2: ef beta omega eps ! If this is susceptibility calculation
  !line 2: emin emax ne eps  ! If this is band calculation
  !line 3: nkx nky nkz     ! If spectral calculation, no use
  !line 4-6 : ACELL
  !line 7: mode            ! Calculation mode
  !  mode=0 : single point calculation
  !  mode=1 : line-mode (band k-path)
  !  mode=2 : plane-mode
  !  mode=3 : whole BZ
  !  Above are diagonal trace only
  !  mode+10 : Full matrix
  !  mode+20 : with self-energy
  !  mode+30 : RPA
  !line 5.. N : depend on mode%10
  !*************************************
  !
  use constants, only : dp, eps4, fin, fout
  use para,      only : inode, para_sync_int, para_sync_real
  use linalgwrap, only : invmat
  !
  implicit none
  !
  character(*) codename
  integer nq_per_sec, iqsec, iq
  integer iqx, iqy, iqz
  real(dp), allocatable :: qbnd_vec(:, :)
  real(dp), dimension(8) :: tt_real
  integer, dimension(2) :: tt_int
  !character dir
  !
  if (inode.eq.0) then
    open(unit=fin, file=trim(codename)//".inp")
    !
    read(fin, *) seed
    !
    read(fin, *) tt_real(1:4)
    read(fin, *) tt_real(5:7)
    do iq=1, 3
      read(fin, *) acell(:, iq)
    enddo
    read(fin, *) tt_real(8)
  endif
  !
  call para_sync_real(tt_real, 8)
  !
  if (codename=='wannchi') then
    mu=tt_real(1)
    beta=tt_real(2)
    omega=tt_real(3)
  else
    nen=nint(tt_real(3))
    allocate(emesh(nen))
    do iq=1, nen
      emesh(iq)=tt_real(1)+(tt_real(2)-tt_real(1))*(iq-1)/(nen-1)
    enddo
  endif
  !
  call para_sync_real(acell, 9)
  bcell=acell
  call invmat(bcell, 3)
  !
  eps=tt_real(4)
  nkx=nint(tt_real(5))
  nky=nint(tt_real(6))
  nkz=nint(tt_real(7))
  mode=nint(tt_real(8))
  level=mode/10
  mode=mod(mode,10)
  !
  if (beta<0) beta=1.d7
  if (eps>eps4.or.eps<eps9) eps=eps4
  !
  select case (mode)
    case (0)
      nqpt=1
      allocate(qvec(1:3, 1:1))
      if (inode.eq.0) read(fin, *) qvec(:, 1)
      call para_sync_real(qvec, 3*nqpt)
    case (1)
      if (inode.eq.0) read(fin, *) tt_int(1:2)
      call para_sync_int(tt_int, 2)
      nqsec=tt_int(1)
      nq_per_sec=tt_int(2)
      nqpt=(nqsec-1)*nq_per_sec+1
      allocate(qbnd_vec(1:3, 1:nqsec))
      allocate(qvec(1:3, nqpt))
      allocate(xq(nqsec))
      allocate(xlabel(nqsec))
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
      call para_sync_real(qvec, 3*nqpt)
    case (2)
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
 END SUBROUTINE
  !
 SUBROUTINE finalize_input
  !
  implicit none
  !
  if (allocated(qvec)) deallocate(qvec)
  if (allocated(emesh)) deallocate(emesh)
  !
 END SUBROUTINE
  !
END MODULE
