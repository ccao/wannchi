SUBROUTINE output_chi(qv, chi)
  !
  use constants, only : dp, fout
  use para,      only : inode
  !
  implicit none
  !
  real(dp) qv(1:3)
  complex(dp) chi
  !
  if (inode.eq.0) then
    write(fout, '(3F12.8,2F22.16)') qv, chi
  endif
  !
END SUBROUTINE
