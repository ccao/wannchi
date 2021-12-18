SUBROUTINE calc_hk(hk, kvec)
  !
  use constants,     only: dp, cmplx_i, twopi, cmplx_0
  use wanndata,      only: nrpt, rvec, norb, weight, ham
  !
  implicit none
  !
  complex(dp), dimension(norb, norb) :: hk
  real(dp), dimension(3) :: kvec
  integer ir
  real(dp)    :: rdotk
  complex(dp) :: fact
  !
  hk(:,:)=cmplx_0
  do ir=1, nrpt
    rdotk=sum(kvec(:)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
    hk=hk+fact*ham(:, :, ir)
  enddo
  !
END SUBROUTINE

SUBROUTINE ham_fix_static()
  !
  use constants,       only: dp
  use wanndata,        only: ham, r000
  use impurity,        only: fulldim, basis_map, restore_lattice, sinf, ncol
  !
  implicit none
  !
  complex(dp), dimension(fulldim, fulldim) :: siglat
  complex(dp), dimension(ncol)             :: soo
  !
  integer ii, jj
  do ii=1, ncol
    soo(ii)=sinf(ii)
  enddo
  !
  call restore_lattice(siglat, soo)
  !
  do ii=1, fulldim
    do jj=1, fulldim
      ham(basis_map(ii), basis_map(jj), r000)=ham(basis_map(ii), basis_map(jj), r000)+siglat(ii, jj)
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE calc_g0(gf, hk, w, ndim, inv)
  !
  use constants,       only: dp
  use linalgwrap,      only: invmat
  !
  implicit none
  !
  integer ndim
  complex(dp), dimension(ndim, ndim) :: gf, hk
  complex(dp) :: w
  logical :: inv ! if true, gf contains G^-1 instead of G
  !
  integer ii
  !
  gf(:,:)=-hk(:,:)
  do ii=1, ndim
    gf(ii, ii)=w+gf(ii, ii)
  enddo
  !
  if (.not.inv) then
    call invmat(gf, ndim)
  endif
  !
END SUBROUTINE

SUBROUTINE calc_corr_realgf(gf, hk, w, inv)
  !
  use constants,        only: dp
  use wanndata,         only: norb
  use impurity,         only: fulldim, basis_map, ncol, restore_lattice, ismatsubara, interpolate_sigma
  use linalgwrap,       only: invmat
  !
  implicit none
  !
  complex(dp), dimension(norb, norb) :: gf, hk
  complex(dp) :: w
  logical :: inv   ! if true, gf contains G^-1 instead of G
  !
  integer ii, jj
  complex(dp), dimension(ncol) :: sigpack
  complex(dp), dimension(fulldim, fulldim) :: sigfull
  !
  if (ismatsubara) then
    write(*, *) "!!! FATAL: self energy is matsubara!"
    stop
  endif
  !
  call interpolate_sigma(sigpack, real(w))
  !
  gf(:,:)=-hk(:,:)
  call restore_lattice(sigfull, sigpack) ! Restore Full self energy
  !
  do ii=1, fulldim
    do jj=1, fulldim
      gf(basis_map(ii), basis_map(jj))=gf(basis_map(ii), basis_map(jj))-sigfull(ii, jj)
    enddo
  enddo
  !
  do ii=1, norb
    gf(ii, ii)=w+gf(ii, ii)
  enddo
  !
  if (.not. inv) then
    call invmat(gf, norb)
  endif
  !
END SUBROUTINE

SUBROUTINE calc_corr_matsgf(gf, hk, iom, inv)
  !
  use constants,        only: dp, cmplx_0, cmplx_i, twopi
  use wanndata,         only: norb
  use input,            only: beta
  use impurity,         only: fulldim, basis_map, nfreq, ncol, restore_lattice, sigma, ismatsubara
  use linalgwrap,       only: invmat
  !
  implicit none
  !
  complex(dp), dimension(norb, norb) :: gf, hk
  integer :: iom
  logical :: inv   ! if true, gf contains G^-1 instead of G
  !
  integer ii, jj
  complex(dp), dimension(ncol) :: sigpack
  complex(dp), dimension(fulldim, fulldim) :: sigfull
  complex(dp) :: w
  !
  if (.not.ismatsubara) then
    write(*, *) "!!! FATAL: self energy is not matsubara!"
    stop
  endif
  !
  w=(iom-0.5d0)*twopi/beta*cmplx_i
  !
  if (iom>nfreq .or. iom<-nfreq+1) then
    sigpack=cmplx_0  ! DO WE NEED TO CORRECT THIS???!!!
  else if (iom<=0) then
    sigpack=conjg(sigma(:,1-iom))
  else
    sigpack=sigma(:, iom)
  endif
  !
  gf(:,:)=-hk(:,:)
  call restore_lattice(sigfull, sigpack) ! Restore Full self energy
  do ii=1, fulldim
    do jj=1, fulldim
      gf(basis_map(ii), basis_map(jj))=gf(basis_map(ii), basis_map(jj))-sigfull(ii, jj)
    enddo
  enddo
  !
  do ii=1, norb
    gf(ii, ii)=w+gf(ii, ii)
  enddo
  !
  if (.not.inv) then
    call invmat(gf, norb)
  endif
  !
END SUBROUTINE

SUBROUTINE calc_partial_hk(hk, kvec, ndim, idx)
  !
  use constants,     only: dp, cmplx_i, twopi, cmplx_0
  use wanndata,      only: nrpt, rvec, norb, weight, ham
  !
  implicit none
  !
  integer ndim
  integer, dimension(ndim) :: idx
  complex(dp), dimension(ndim, ndim) :: hk
  real(dp), dimension(3) :: kvec
  !
  integer ir, ii, jj
  real(dp)    :: rdotk
  complex(dp) :: fact
  !
  complex(dp), dimension(norb, norb) :: htmp
  !
  htmp(:,:)=cmplx_0
  do ir=1, nrpt
    rdotk=sum(kvec(:)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
    htmp=htmp+fact*ham(:, :, ir)
  enddo
  !
  do ii=1, ndim
    do jj=1, ndim
      hk(ii, jj)=htmp(idx(ii), idx(jj))
    enddo
  enddo
  !hk(:,:)=htmp(idx(:),idx(:))
  !
END SUBROUTINE
