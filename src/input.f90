MODULE input
  !
  use constants
  !
  implicit none
  !
  integer nqpt
  real(dp), allocatable :: qvec(:, :)
  real(dp) temp, omega, eps
  character(len=80) seed
  !
CONTAINS
  !
 SUBROUTINE read_input
  !
  !************ INPUT FILE *************
  !** file name: wannchi.inp
  !line 1: seed name
  !line 2: ef T omega epsilon
  !line 3: nqx nqy nqz
  !line 4: mode
  !line 5: ...
  !*************************************
  !
  use constants, only : dp, eps4, fin
  use para,      only : inode, para_sync_int, para_sync_real
  use banddata,  only : nkpt, nkx, nky, nkz, ef
  !
  implicit none
  !
  integer mode
  integer nqsec, nq_per_sec, iqsec, iq
  integer iqx, iqy, iqz
  real(dp), allocatable :: qbnd_vec(:, :)
  real(dp), dimension(4) :: tt_real
  integer, dimension(4)  :: tt_int
  !character dir
  !
  if (inode.eq.0) then
    open(unit=fin, file="wannchi.inp")
    !
    read(fin, *) seed
    !
    read(fin, *) tt_real(1:4)
    read(fin, *) tt_int(1:3)
    read(fin, *) mode
  endif
  !
  call para_sync_int(tt_int, 4)
  call para_sync_real(tt_real, 4)
  !
  ef=tt_real(1)
  temp=tt_real(2)
  omega=tt_real(3)
  eps=tt_real(4)

  nkx=tt_int(1)
  nky=tt_int(2)
  nkz=tt_int(3)
  mode=tt_int(4)
  !
  nkpt=nkx*nky*nkz
  !
  if (temp<0) temp=0.d0
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
      if (inode.eq.0) then
        do iqsec=1, nqsec
          read(fin, *) qbnd_vec(:, iqsec)
        enddo
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
  deallocate(qvec)
  !
 END SUBROUTINE
  !
END MODULE
