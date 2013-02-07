SUBROUTINE read_input
  !
  use banddata
  !
  implicit none
  !
  open(unit=fin, file="input")
  read(fin, *) ef
  read(fin, *) nkx, nky, nkz
  close(unit=fin)
  !
  nkpt=nkx*nky*nkz
  !
  allocate(eig(1:nkpt, 1:nbnd))
  allocate(egv(1:nkpt, 1:nbnd, 1:nbnd))
  !
END SUBROUTINE
