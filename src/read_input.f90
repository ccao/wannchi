SUBROUTINE read_input
  !
  use para
  use banddata
  !
  implicit none
  !
  if (inode.eq.0) then
    open(unit=fin, file="input")
    read(fin, *) ef
    read(fin, *) nkx, nky, nkz
    close(unit=fin)
  endif
  !
  CALL para_sync(nbnd)
  CALL para_sync(ef)
  CALL para_sync(nkx)
  CALL para_sync(nky)
  CALL para_sync(nkz)
  !
  nkpt=nkx*nky*nkz
  !
  allocate(eig(1:nbnd, 1:nkpt))
  allocate(egv(1:nbnd, 1:nbnd, 1:nkpt))
  !
END SUBROUTINE
