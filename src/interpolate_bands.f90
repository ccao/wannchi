include 'lapack.f90'

SUBROUTINE interpolate_bands
  !
  use para
  use lapack95,  only : heev
  use constants, only : dp, twopi, cmplx_0, cmplx_i, stdout
  use wanndata,  only : rvec, ham, weight, nrpt, norb
<<<<<<< HEAD
  use banddata,  only : nkpt, nkx, nky, nkz, egv, eig
  use input,     only : code
=======
  use banddata,  only : nkpt, nkx, nky, nkz, xi, eig, occ, ef, kvec
>>>>>>> New modulized version
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
  integer first_k, last_k
  !
  allocate(work(1:norb, 1:norb))
  allocate(e(1:norb))
  !
  if(inode.eq.0) then
    write(stdout, *) " # Starting interpolation of the original states to"
    write(stdout, *) " # ", nkx, "x", nky, "x", nkz, " K-mesh"
  endif
  !
  first_k=inode*nkpt/nnode+1
  last_k=(inode+1)*nkpt/nnode
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
<<<<<<< HEAD
    if (code.eq.0) &
     & egv(:, :, ik)=work(:,:)
=======
>>>>>>> New modulized version
    eig(:, ik)=e(:)
    !
    do ii=1, norb
      occ(ii, ik)=calc_occ(eig(ii, ik))
      do io=1, norb
        do jo=1, norb
          xi(io, jo, ii, ik)=work(io, ii)*conjg(work(jo, ii))
        enddo
      enddo
    enddo
    !
  enddo ! ik
  !
<<<<<<< HEAD
  if (code.eq.0) &
     & CALL para_merge(egv, norb, norb, nkpt)
=======
  CALL para_merge(xi, norb, norb, norb, nkpt)
>>>>>>> New modulized version
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
