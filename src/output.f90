SUBROUTINE output_header(fd)
  !
  use constants, only : dp
  use input,     only : emesh, nen, nqpt, qvec, mode, nqsec, xq, xlabel
  implicit none
  !
  integer fd
  integer ii, jj
  integer nq_per_sec
  real(dp), dimension(nqpt) :: xxq
  if (mode==0) then
    write(fd, '(A,3F9.4)') "# Spectral function for kpt: ", qvec(:, 1)
  else if (mode==1) then
    nq_per_sec=(nqpt-1)/(nqsec-1)
    do ii=1, nqsec
      write(fd, '(1F14.9,15X,A)') xq(ii), xlabel(ii)
    enddo
    !
    write(fd, '(2I10)') nqpt, nen
    !
    do ii=1, nqsec-1
      do jj=1, nq_per_sec
        xxq((ii-1)*nq_per_sec+jj)=xq(ii)+(xq(ii+1)-xq(ii))*(jj-1)/nq_per_sec
      enddo
    enddo
    xxq(nqpt)=xq(nqsec)
    write(fd, '(10F14.9)') xxq(:)
    write(fd, '(10F14.9)') emesh(:)
  endif
  !
END SUBROUTINE

SUBROUTINE output_spectral(spec, fd, nq)
  !
  use constants, only : dp
  use input,     only : emesh, nen
  implicit none
  !
  integer fd, nq
  real(dp), dimension(nen) :: spec
  !
  integer ii
  !
  if (nq.eq.1) then
    do ii=1, nen
      write(fd, '(2E18.10)') emesh(ii), spec(ii)
    enddo
  else
    write(fd, '(10E18.10)') spec(:)
  endif
  !
END SUBROUTINE

SUBROUTINE output_chi(chi, fd, nq)
  !
  use constants, only : dp
  use input,     only : emesh, nen
  implicit none
  !
  integer fd, nq
  complex(dp), dimension(nen) :: chi
  !
  real(dp), dimension(nen) :: spec
  integer ii
  !
  if (nq.eq.1) then
    do ii=1, nen
      write(fd, '(3E18.10)') emesh(ii), chi(ii)
    enddo
  else
    spec(:)=real(chi(:))
    write(fd, '(10E18.10)') spec(:)
  endif
  !
END SUBROUTINE
