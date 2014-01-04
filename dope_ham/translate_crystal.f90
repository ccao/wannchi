LOGICAL FUNCTION isinteger(a)
  !
  use constants
  !
  IMPLICIT NONE
  !
  real(dp) a
  !
  if ( abs(a-nint(a)).le.1e-5 ) then
    isinteger=.true.
  else
    isinteger=.false.
  endif
  return
  !
END FUNCTION

SUBROUTINE translate_crystal(trans_vec)
  !
  use constants
  use orbitals
  use dopedata
  !
  IMPLICIT NONE
  !
  real(dp) :: trans_vec(1:3)
  real(dp) :: pos(1:3)
  real(dp) :: dx(1:3)
  integer ii, jj
  logical isinteger, found
  !
  do ii=1, nat
    found=.false.
    pos(:)=tau(ii, :)+trans_vec(:)
    do jj=1, nat
      dx(:)=pos(:)-tau(jj,:)
      if (isinteger(dx(1)).and.isinteger(dx(2)).and.isinteger(dx(3))) then
        if(at_sp(ii).ne.at_sp(jj)) then
          write(stdout, *) " !!! ERROR: Internal error # 1"
          stop
        endif
        new_to_old(jj)=ii
        new_to_old_dR(jj,:)=nint(dx(:))
        found=.true.
        exit ! loop jj
      endif
    enddo
    if (.not.found) then
      write(stdout, *) " !!! ERROR: Cannot identify translated atom:"
      write(stdout, *) " translated coordinate:"
      write(stdout, *) pos(:)
      write(stdout, *) " Atomic specie:", at_sp(ii)
      stop
    endif
  enddo
  !
END SUBROUTINE
