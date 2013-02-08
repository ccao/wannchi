include 'lapack.f90'

SUBROUTINE interpolate_bands
  !
  use para
  use lapack95,  only : heev
  use constants, only : dp, twopi, cmplx_0, cmplx_i, stdout
  use wanndata,  only : rvec, ham, weight, nrpt, norb
  use banddata,  only : nkpt, nkx, nky, nkz, egv, eig
  !
  implicit none
  !
  real(dp) kvec(1:3)
  real(dp) rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:,:)
  real(dp), allocatable :: e(:)
  !
  integer ir, ik, info
  integer first_k, last_k
  !
  allocate(work(1:norb, 1:norb))
  allocate(e(1:norb))
  !
  if(inode.eq.0) then
    write(stdout, *) "Starting interpolation of the original states to"
    write(stdout, *) nkx, "x", nky, "x", nkz, " K-mesh"
  endif
  !
  first_k=inode*nkpt/nnode+1
  last_k=(inode+1)*nkpt/nnode
  do ik=first_k, last_k
    CALL map_to_kpt(kvec, ik, nkx, nky, nkz)
    !
    work(:,:)=cmplx_0
    do ir=1, nrpt
      rdotk=SUM(kvec(:)*rvec(:, ir))
      fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
      work(:,:)=work(:,:)+ham(:,:,ir)*fact
    enddo ! ir
    !
    call heev(work, e, 'V', 'U', info)
    !
    egv(:, :, ik)=work(:,:)
    eig(:, ik)=e(:)
  enddo ! ik
  !
  CALL para_merge(egv, norb, norb, nkpt)
  CALL para_merge(eig, norb, nkpt)
  !
  if(inode.eq.0) then
    write(stdout, *) "Done..."
  endif
  !
  deallocate(work)
  !
END SUBROUTINE
