include 'lapack.f90'

SUBROUTINE interpolate_bands
  !
  use para,      only : inode, nnode, para_merge, first_k, last_k
  use lapack95,  only : heev
  use constants, only : dp, twopi, cmplx_0, cmplx_i, stdout
  use wanndata,  only : rvec, ham, weight, nrpt, norb
  use banddata,  only : nkpt, nkx, nky, nkz, egv, eig, occ, ef, kvec
  !
  implicit none
  !
  real(dp) rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:,:)
  real(dp), allocatable :: e(:)
  real(dp) calc_occ
  !
  integer ir, ik, info, io, jo, ii
  !
  allocate(work(1:norb, 1:norb))
  allocate(e(1:norb))
  !
  if(inode.eq.0) then
    write(stdout, *) " # Starting interpolation of the original states to"
    write(stdout, *) " # ", nkx, "x", nky, "x", nkz, " K-mesh"
    write(stdout, *) " # Fermi level :", ef
  endif
  !
  do ik=first_k, last_k
    work(:,:)=cmplx_0
    do ir=1, nrpt
      rdotk=SUM(kvec(:, ik)*rvec(:, ir))
      fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
      work(:,:)=work(:,:)+ham(:,:,ir)*fact
    enddo ! ir
    !
    call heev(work, e, 'V', 'U', info)
    !
    eig(:, ik)=e(:)
    egv(:, :, ik)=work(:, :)  ! egv(ii, io, ik): ii\th element eigenvector of eigenvalue (io, ik)
    do ii=1, norb
      occ(ii, ik)=calc_occ(eig(ii, ik)-ef)
    enddo
    !
  enddo ! ik
  !
  CALL para_merge(occ, norb, nkpt)
  CALL para_merge(egv, norb, norb, nkpt)
  CALL para_merge(eig, norb, nkpt)
  !
  if(inode.eq.0) then
    write(stdout, *) " # Done..."
  endif
  !
  deallocate(work)
  deallocate(e)
  !
END SUBROUTINE
