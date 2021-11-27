!
!   impurity.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2021 CC. All rights reserved.
!
MODULE impurity
  !
  use constants
  !
  implicit none
  !
  logical ismatsubara
  ! the self energy is matsubara frequency
  integer nimp
  ! number of impurities
  integer, dimension(:), allocatable :: ndim
  ! Dimension of each impurity
  integer :: fulldim
  ! Dimension of the full impurity problem (=sum(ndim))
  integer, dimension(:), allocatable :: nnlow
  ! Lower bound of ii-th impurity index is nnlow(ii)
  integer, dimension(:, :, :), allocatable :: sigind
  ! mapping from sig.inp to actual self energy matrix
  complex(dp), dimension(:, :, :), allocatable :: Utrans
  ! Transform matrix from DMFT basis to spherical harmonix
  complex(dp), dimension(:, :), allocatable :: Ufull
  ! Full transform matrix from DMFT basis to lattice basis
  integer, dimension(:), allocatable :: basis_map
  ! Mapping from impurity to lattice
  ! i-th orbital in impurity problem is basis_map(i)-th orbital in lattice problem
  !
  integer ncol
  integer nfreq
  ! self energy has ncolumns and nfrequencies 
  ! Sigma(ncol, nfreq)
  real(dp), dimension(:), allocatable :: sinf
  ! soo-vdc sinf(ncol)
  complex(dp), dimension(:), allocatable :: omega
  ! frequency mesh
  complex(dp), dimension(:, :), allocatable :: sigma
  ! Actual self energy
  !
  contains
  !
  subroutine restore_lattice(siglat, sigpack)
    !
    use blas95,     only : gemm
    implicit NONE
    !
    complex(dp), dimension(fulldim, fulldim) :: siglat
    complex(dp), dimension(ncol) :: sigpack
    !
    complex(dp), dimension(fulldim, fulldim) :: sigtmp
    !
    integer ii, io1, io2
    siglat=cmplx_0
    do ii=1, nimp
      do io1=1, ndim(ii)
        do io2=1, ndim(ii)
          if (sigind(io1, io2, ii)>0) siglat(nnlow(ii)+io1, nnlow(ii)+io2, ii)=sigpack(sigind(io1, io2, ii))
        enddo
      enddo
    enddo
    !
    call gemm(siglat, Utrans, sigtmp, 'N', 'N')
    call gemm(Utrans, sigtmp, siglat, 'T', 'N')
    !
  end subroutine

  subroutine finalize_impurity()
    !
    if (allocated(ndim)) deallocate(ndim)
    if (allocated(nnlow)) deallocate(nnlow)
    if (allocated(sigind)) deallocate(sigind)
    if (allocated(Utrans)) deallocate(Utrans)
    if (allocated(basis_map)) deallocate(basis_map)
    if (allocated(sinf)) deallocate(sinf)
    if (allocated(omega)) deallocate(omega)
    if (allocated(sigma)) deallocate(sigma)
    !
  end subroutine

  subroutine init_impurity()
    !
    use para,    only: inode, para_sync
    implicit none
    !
    character(len=80) :: fn_indmfl, fn_siginp
    !
    if (inode.eq.0) then
      open(unit=fin, file='impurity.inp')
      read(fin, *) ncol, nfreq
      read(fin, *) fn_indmfl
      read(fin, *) fn_siginp
      close(unit=fin)
    endif
    !
    call para_sync(ncol)
    call para_sync(nfreq)
    !
    allocate(omega(nfreq))
    allocate(sigma(ncol, nfreq))
    allocate(sinf(ncol))
    !
    call read_indmfl(fn_indmfl)
    call read_siginp(fn_siginp)
    !
  end subroutine
  !
  subroutine read_siginp(fname)
    !
    use para,    only: inode, para_sync
    implicit none
    !
    character(len=80) :: fname
    !
    real(dp), dimension(:), allocatable :: aa
    !
    integer ii, jj
    !
    if (inode.eq.0) then
      !
      allocate(aa(2*ncol+1))
      !
      open(unit=fin, file=trim(fname))
      read(fin, *) sinf
      read(fin, *) aa(1:ncol)
      sinf(:)=sinf(:)-aa(1:ncol)
      !
      do ii=1, nfreq
        read(fin, *) aa
        if (ismatsubara) then
          omega(ii)=aa(1)*cmplx_i
        else
          omega(ii)=aa(1)
        endif
        do jj=1, ncol
          sigma(jj, ii)=aa(2*jj)+aa(2*jj+1)*cmplx_i
        enddo
      enddo
      !
      close(unit(fin))
      !
      deallocate(aa)
      !
    endif
    !
    call para_sync(omega, nfreq)
    call para_sync(sinf, ncol)
    call para_sync(sigma, ncol*nfreq)
    !
  end subroutine
  !
  subroutine read_indmfl(fname)
    !
    use para,    only: inode, para_sync
    implicit none
    !
    character(len=80) :: fname
    integer t1, t2, t3
    integer maxdim
    integer ii, jj
    real(dp), dimension(:), allocatable :: aa
    !
    if (inode.eq.0) then
      open(unit=fin, file=trim(fname))
      read(fin, *)     !  # hybridization band index nemin and nemax, renormalize for interstitials, projection type
      read(fin, *) t1  !  # matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)
    endif
    !
    call para_sync(t1)
    ismatsubara=(t1.eq.1)
    !
    if (inode.eq.0) then
      read(fin, *) nimp ! # number of correlated atoms
      do ii=1, nimp
        read(fin, *)   !  # iatom, nL, locrot
        read(fin, *)   !  # L, qsplit, cix
      enddo
      !
      read(fin, *)     !  #================ # Siginds and crystal-field transformations for correlated orbitals ================
      read(fin, *) nimp, maxdim, t1
      !
    endif
    !
    call para_sync(nimp)
    call para_sync(maxdim)
    allocate(ndim(nimp))
    allocate(nnlow(nimp))
    allocate(sigind(maxdim, maxdim, nimp))
    allocate(Utrans(maxdim, maxdim, nimp))
    !
    sigind(:, :, :)=-1
    Utrans(:, :, :)=cmplx_0
    !
    if (inode.eq.0) then
      !
      allocate(aa(2*maxdim))
      !                   # Number of independent kcix blocks, max dimension, max num-independent-components
      do ii=1, nimp
        read(fin, *) t1, t2, t3
        !                 # cix-num, dimension, num-independent-components
        ndim(ii)=t2
        read(fin, *)   !  #---------------- # Independent components are --------------
        read(fin, *)   !  '5/2' '7/2'
        read(fin, *)   !  #---------------- # Sigind follows --------------------------
        do jj=1, ndim(ii)
          read(fin, *) sigind(1:ndim(ii), jj, ii)
        enddo
        !
        read(fin, *)   !  #---------------- # Transformation matrix follows -----------
        do jj=1, ndim(ii)
          read(fin, *) aa(1:2*ndim(ii))
          do kk=1, ndim(ii)
            Utrans(kk, jj, ii)=aa(2*kk-1)+aa(2*kk)*cmplx_i
          enddo
        enddo
      enddo
      ！
      close(unit=fin)
    endif
    !
    call para_sync(ndim, nimp)
    call para_sync(sigind, maxdim*maxdim*nimp)
    call para_sync(Utrans, maxdim*maxdim*nimp)
    !
    fulldim=sum(ndim)
    nnlow(1)=1
    do ii=2, nimp
      nnlow(ii)=nnlow(ii-1)+ndim(ii-1)
    enddo
    !
  end subroutine
  !
END MODULE