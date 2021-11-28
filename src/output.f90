SUBROUTINE output_header(fd)
  !
  use constants, only : dp
  use input,     only : emesh, nen, nqpt, qvec, mode
  implicit none
  !
  integer fd
  integer ii
  real(dp), dimension(nqpt) :: xq
  if (mode==0) then
    write(fd, '(A,3F9.4)') "# Spectral function for kpt: ", qvec(:, 1)
  else if (mode==1) then
    do ii=1, nqpt
      xq(ii)=(ii-1.0)*1.d0/(nqpt-1.0)
    enddo
    !
    write(fd, '(1F9.4,15X,A)') xq(1), 'G'
    write(fd, '(1F9.4,15X,A)') xq(nqpt), 'Z'
    !
    write(fd, '(2I10)') nqpt, nen
    !
    write(fd, '(10F14.9)') xq(:)
    write(fd, '(10F14.9)') emesh(:)
  endif
  !
END SUBROUTINE

SUBROUTINE output_spectral(spec, ene, nen, fd, nq)
  !
  use constants, only : dp
  implicit none
  !
  integer nen, fd, nq
  real(dp), dimension(nen) :: spec, ene
  !
  integer ii
  !
  if (nq.eq.1) then
    do ii=1, nen
      write(fd, '(2F18.12)') ene(ii), spec(ii)
    enddo
  else
    write(fd, '(10F14.9)') spec(:)
  endif
  !
END SUBROUTINE