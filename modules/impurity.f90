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
          if (sigind(io1, io2, ii)>0) siglat(nnlow(ii)+io1, nnlow(ii)+io2)=sigpack(sigind(io1, io2, ii))
        enddo
      enddo
    enddo
    !
    call zgemm('N', 'N', fulldim, fulldim, fulldim, cmplx_1, siglat, fulldim, Utrans, fulldim, cmplx_0, sigtmp, fulldim)
    call zgemm('T', 'N', fulldim, fulldim, fulldim, cmplx_1, Utrans, fulldim, sigtmp, fulldim, cmplx_0, siglat, fulldim)
    !call gemm(siglat, Utrans, sigtmp, 'N', 'N')
    !call gemm(Utrans, sigtmp, siglat, 'T', 'N')
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
    use para,    only: inode, para_sync_int
    implicit none
    !
    character(len=80) :: fn_indmfl, fn_siginp
    integer, dimension(2) :: tt
    !
    if (inode.eq.0) then
      open(unit=fin, file='impurity.inp')
      read(fin, *) tt
      read(fin, *) fn_indmfl
      read(fin, *) fn_siginp
      close(unit=fin)
    endif
    !
    call para_sync_int(tt, 2)
    ncol=tt(1)
    nfreq=tt(2)
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
    use para,    only: inode, para_sync_int, para_sync_real, para_sync_cmplx
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
      write(stdout, *) ' Reading Self Energies from '//trim(fname)
      open(unit=fin, file=trim(fname))
      read(fin, *) sinf
      read(fin, *) aa(1:ncol)
      sinf(:)=sinf(:)-aa(1:ncol)
      write(stdout, *) '   Static part (Sigma(infty)-Vdc) are:'
      do ii=1, ncol
        write(stdout, '(1F14.9,2X)', advance='no') sinf(ii)
      enddo
      write(stdout, *)
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
      close(fin)
      !
      deallocate(aa)
      !
    endif
    !
    call para_sync_cmplx(omega, nfreq)
    call para_sync_real(sinf, ncol)
    call para_sync_cmplx(sigma, ncol*nfreq)
    !
  end subroutine
  !
  subroutine read_indmfl(fname)
    !
    use para,    only: inode, para_sync_int, para_sync_cmplx
    implicit none
    !
    character(len=80) :: fname
    integer, dimension(3) :: tt
    integer maxdim
    integer ii, jj, kk
    real(dp), dimension(:), allocatable :: aa
    !
    if (inode.eq.0) then
      open(unit=fin, file=trim(fname))
      write(stdout, *) ' Reading impurity definitions from'//trim(fname)
      read(fin, *)     !  # hybridization band index nemin and nemax, renormalize for interstitials, projection type
      read(fin, *) tt(1) ! # matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)
      read(fin, *) tt(2) ! # number of correlated atoms
      do ii=1, tt(2)
        read(fin, *)   !  # iatom, nL, locrot
        read(fin, *)   !  # L, qsplit, cix
      enddo
      !
      read(fin, *)     !  #================ # Siginds and crystal-field transformations for correlated orbitals ================
      read(fin, *) tt(2:3)
      !
    endif
    !
    call para_sync_int(tt, 3)
    ismatsubara=(tt(1).eq.1)
    nimp=tt(2)
    maxdim=tt(3)
    !
    allocate(ndim(nimp))
    allocate(nnlow(nimp))
    allocate(sigind(maxdim, maxdim, nimp))
    allocate(Utrans(maxdim, maxdim, nimp))
    !
    sigind(:, :, :)=-1
    Utrans(:, :, :)=cmplx_0
    !
    if (inode.eq.0) then
      if (ismatsubara) then
        write(stdout, *) '   Self energy is Matsubara'
      else
        write(stdout, *) '   Self energy is real-frequency'
      endif
      write(stdout, '(A,1I3,A)') '   There are ', nimp, ' impurities'
      !write(stdout, '(A,1I3)')   '     Maximum dimension:', maxdim
      !
      allocate(aa(2*maxdim))
      !                   # Number of independent kcix blocks, max dimension, max num-independent-components
      write(stdout, '(A)')       '     each with dimensions:'
      do ii=1, nimp
        read(fin, *) tt
        !                 # cix-num, dimension, num-independent-components
        ndim(ii)=tt(2)
        write(stdout, '(1I4)', advance='no') ndim(ii)
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
      !
      write(stdout, *)
      !
      close(unit=fin)
    endif
    !
    call para_sync_int(ndim, nimp)
    call para_sync_int(sigind, maxdim*maxdim*nimp)
    call para_sync_cmplx(Utrans, maxdim*maxdim*nimp)
    !
    fulldim=sum(ndim)
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I3)')   '   Full impurity dimension:', fulldim
    endif
    nnlow(1)=1
    do ii=2, nimp
      nnlow(ii)=nnlow(ii-1)+ndim(ii-1)
    enddo
    !
  end subroutine
  !
END MODULE