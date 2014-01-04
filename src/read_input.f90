SUBROUTINE read_input
  !
  use para
  use banddata
  use input
  !
  implicit none
  !
  integer iseg
  !
  if (inode.eq.0) then
    open(unit=fin, file="input")
    read(fin, *) ef
    read(fin, *) nkx, nky, nkz
    if (code.eq.0) then
      read(fin, *) mode
      read(fin, *) temp, omega
      if (mode.eq.1) then
        read(fin, *) nqseg, nqbnd
      else
        read(fin, *) hubbard_u, hubbard_j, hubbard_v, hubbard_jp
      endif
    else  ! code=1
      read(fin, *) sigma
    endif
  endif
  !
  CALL para_sync(nbnd)
  CALL para_sync(ef)
  CALL para_sync(nkx)
  CALL para_sync(nky)
  CALL para_sync(nkz)
  CALL para_sync(mode)
  CALL para_sync(sigma)
  CALL para_sync(nqseg)
  CALL para_sync(nqbnd)
  CALL para_sync(temp)
  CALL para_sync(omega)
  CALL para_sync(hubbard_u)
  CALL para_sync(hubbard_j)
  CALL para_sync(hubbard_v)
  CALL para_sync(hubbard_jp)
  !
  if ( abs(hubbard_u)+abs(hubbard_j)+abs(hubbard_v)+abs(hubbard_jp) > 0.d0 ) then
    lrpa=.true.
    if (inode.eq.0) write(*, *) "# !!!!  None-zero Hubbard Parameters, will do RPA!!!"
  endif
  !
  if(mode.eq.1) then
    allocate(bnd_q(1:nqseg+1, 1:3))
    if(inode.eq.0) then
      do iseg=1, nqseg+1
        read(fin, *) bnd_q(iseg, :)
      enddo
    endif
  CALL para_sync(bnd_q, nqseg+1, 3)
  endif
  !
  nkpt=nkx*nky*nkz
  !
  allocate(eig(1:nbnd, 1:nkpt))
  if (code.eq.0)  &
     & allocate(egv(1:nbnd, 1:nbnd, 1:nkpt))
  !
  close(unit=fin)
  !
END SUBROUTINE
