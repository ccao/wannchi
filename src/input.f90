MODULE input
  !
  use constants
  !
  implicit none
  !
<<<<<<< HEAD
  integer code   ! 0: wannchi 1: wanndos
  integer mode   ! Mode for calculation, 0 (default): plane; 1: band
  real(dp), allocatable :: bnd_q(:, :) ! special Q-points
  real(dp) temp  ! Temperature
  real(dp) omega ! omega ( real axis energy )
  real(dp) hubbard_u, hubbard_j, hubbard_v, hubbard_jp
  integer nqseg   ! Num of q-points per seg
  integer nqbnd   ! Num of special Q points
  !
  logical lrpa
  !
CONTAINS

SUBROUTINE finalize_input()
  !
  implicit none
  !
  if(allocated(bnd_q)) deallocate(bnd_q)
  !
END SUBROUTINE

SUBROUTINE get_qvec(qv, iq)
  !
  implicit none
  !
  integer iq
  real(dp) qv(1:3)
  integer iseg, iqb
  !
  iqb=mod(iq-1,nqbnd)
  iseg=(iq-1-iqb)/nqbnd
  !
  qv(:)=bnd_q(iseg+1,:)+iqb*(bnd_q(iseg+2,:)-bnd_q(iseg+1,:))/nqbnd
  !
END SUBROUTINE

END MODULE

=======
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
  !line 2: fermi level
  !line 3: nqx nqy nqz
  !line 4: temperature
  !line 5: omega
  !line 6: epsilon
  !line 7: mode
  !line 8: ...
  !*************************************
  !
  use constants, only : dp, eps4, fin
  use para
  use banddata,  only : nkpt, nkx, nky, nkz, ef
  !
  implicit none
  !
  integer mode
  integer nqsec, nq_per_sec, iqsec, iq
  integer iqx, iqy, iqz
  real(dp), allocatable :: qbnd_vec(:, :)
  character dir
  !
  if (inode.eq.0) then
    open(unit=fin, file="wannchi.inp")
    !
    read(fin, *) seed
    !
    read(fin, *) ef
    read(fin, *) nkx, nky, nkz
    read(fin, *) temp
    read(fin, *) omega
    read(fin, *) eps
    read(fin, *) mode
    !
    nkpt=nkx*nky*nkz
  endif
  !
  CALL para_sync(ef)
  CALL para_sync(nkx)
  CALL para_sync(nky)
  CALL para_sync(nkz)
  CALL para_sync(temp)
  CALL para_sync(omega)
  CALL para_sync(eps)
  CALL para_sync(mode)
  CALL para_sync(nkpt)
  !
  if (temp<0) temp=0.d0
  if (eps>eps4.or.eps<eps9) eps=eps4
  !
  select case (mode)
    case (0)
      nqpt=1
      allocate(qvec(1:3, 1:1))
      if (inode.eq.0) read(fin, *) qvec(:, 1)
      CALL para_sync(qvec, 3, nqpt)
    case (1)
      read(fin, *) nqsec, nq_per_sec
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
      CALL para_sync(qvec, 3, nqpt)
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
>>>>>>> New modulized version
