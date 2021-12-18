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
  integer, parameter :: max_num_imp=10
  !
  logical ismatsubara
  ! the self energy is matsubara frequency
  integer nimp
  ! number of impurities
  integer, dimension(max_num_imp) :: ndim
  ! Dimension of each impurity (maximally 10 impurities for now)
  integer :: fulldim
  ! Dimension of the full impurity problem (=sum(ndim))
  integer, dimension(:), allocatable :: nnlow
  ! Lower bound of ii-th impurity index is nnlow(ii)
  integer, dimension(:, :, :), allocatable :: sigind
  ! mapping from sig.inp to actual self energy matrix
  complex(dp), dimension(:, :, :), allocatable :: Utrans
  ! Transform matrix from DMFT basis to cubic harmonix
  integer, dimension(:), allocatable :: basis_map
  ! Mapping from impurity to lattice
  ! i-th orbital in impurity problem is basis_map(i)-th orbital in lattice problem
  !
  character(len=80) sigfile
  integer ncol
  integer nfreq
  ! self energy has ncolumns and nfrequencies 
  ! Sigma(ncol, nfreq)
  real(dp), dimension(:), allocatable :: sinf
  ! soo-vdc sinf(ncol)
  real(dp), dimension(:), allocatable :: omega
  ! frequency mesh, used only in real-frequency calculations
  complex(dp), dimension(:, :), allocatable :: sigma
  ! Actual self energy
  namelist /IMPDEF/ ismatsubara, sigfile, ncol, nfreq, nimp, ndim
  !
  contains
  !
  subroutine init_impurity(beta)
    !
    use para,    only: inode, para_sync_int, para_sync_cmplx
    use symmetry_module, only: generate_Ylm2C, symmetry, rotate_J, rotate_cubic, init_symm, inverse_symm
    !
    implicit none
    !
    real(dp)  :: beta
    !
    character(len=80) :: line
    integer ii, jj, kk, maxdim, corr_l
    integer, dimension(4) :: tt
    real(dp), dimension(:, :), allocatable :: locrot
    real(dp), dimension(:), allocatable :: aa
    complex(dp), dimension(:, :), allocatable :: umat, mtmp, locrot_j
    TYPE(symmetry) :: symm
    !
    ! impurity.inp:
    ! &IMPDEF
    !   ismatsubara=.true.
    !   sigfile='sig.inp'
    !   ncol=2
    !   nfreq=40000
    !   nimp=2
    !     ndim(1)=14
    !     ndim(2)=14
    ! /
    ! Impurities
    ! BasisMap
    !   1 2 3 4 5 6 7  61 62 63 64 65 66 67
    ! LocalAxis
    !   1.0 0.0 0.0
    !   0.0 1.0 0.0
    !   0.0 0.0 1.0
    ! Sigind
    !   1 0 0 0 0 0 0 0 0 0 0 0 0 0
    !   0 1 0 0 0 0 0 0 0 0 0 0 0 0
    !   0 0 1 0 0 0 0 0 0 0 0 0 0 0
    !   0 0 0 1 0 0 0 0 0 0 0 0 0 0
    !   0 0 0 0 1 0 0 0 0 0 0 0 0 0
    !   0 0 0 0 0 1 0 0 0 0 0 0 0 0
    !   0 0 0 0 0 0 2 0 0 0 0 0 0 0
    !   0 0 0 0 0 0 0 2 0 0 0 0 0 0
    !   0 0 0 0 0 0 0 0 2 0 0 0 0 0
    !   0 0 0 0 0 0 0 0 0 2 0 0 0 0
    !   0 0 0 0 0 0 0 0 0 0 2 0 0 0
    !   0 0 0 0 0 0 0 0 0 0 0 2 0 0
    !   0 0 0 0 0 0 0 0 0 0 0 0 2 0
    !   0 0 0 0 0 0 0 0 0 0 0 0 0 2
    ! Transformation Matrix
    ! ....
    !
    ! DEFAULTS...
    ismatsubara = .true.
    sigfile='sig.inp'
    ncol=1
    nfreq=1
    nimp=1
    ndim(:)=0
    !
    if (inode.eq.0) then
      open(unit=fin, file='impurity.inp')
      read(nml=IMPDEF, unit=fin)
      tt(1)=ncol
      tt(2)=nfreq
      tt(3)=nimp
      if (ismatsubara) then
        tt(4)=1
      else
        tt(4)=0
      endif
      write(stdout, *) " Reading impurity definitions."
    endif
    !
    call para_sync_int(tt, 4)
    ncol=tt(1)
    nfreq=tt(2)
    nimp=tt(3)
    ismatsubara=(tt(4).eq.1)
    call para_sync_int(ndim, nimp)
    !
    allocate(omega(nfreq))
    allocate(sigma(ncol, nfreq))
    allocate(sinf(ncol))
    allocate(nnlow(nimp))
    !
    maxdim=0
    fulldim=0
    do ii=1, nimp
      if (ndim(ii)>maxdim) maxdim=ndim(ii)
      nnlow(ii)=fulldim+1
      fulldim=fulldim+ndim(ii)
    enddo
    !
    allocate(basis_map(fulldim))
    allocate(Utrans(maxdim, maxdim, nimp))
    allocate(sigind(maxdim, maxdim, nimp))
    basis_map(:)=0
    Utrans(:,:,:)=cmplx_0
    sigind(:,:,:)=0
    !
    if (inode.eq.0) then
      allocate(locrot(3, 3), aa(2*maxdim))
      write(stdout, '(1I5,A)') nimp, ' impurities will be read'
      !
      do ii=1, nimp
        !
        read(fin, *) line
        if (trim(line(1:4))/='IMP') then
          write(*, *) "!!! FATAL ERROR: Incorrect format IMP"
          stop
        endif

        read(fin, *) line
        if (trim(line)/='BasisMap') then
          write(*, *) "!!! FATAL ERROR: Incorrect format BasisMap"
          stop
        endif
        read(fin, *) basis_map(nnlow(ii):nnlow(ii)+ndim(ii)-1)
        write(stdout, '(A,1I4,A,1I4)') "Dimension of impurity ", ii, " :", ndim(ii)
        write(stdout, *) "  Basis are :"
        do jj=1, ndim(ii)
          write(stdout, '(1I5)', advance='no') basis_map(nnlow(ii)+jj-1)
        enddo
        write(stdout, *)
        !
        read(fin, *) line
        if (trim(line)/='LocalAxis') then
          write(*, *) "!!! FATAL ERROR: Incorrect format LocalAxis"
          stop
        endif
        do jj=1, 3
          read(fin, *) locrot(:, jj)
        enddo
        !
        read(fin, *) line
        if (trim(line(1:9))/='Sigind') then
          write(*, *) "!!! FATAL ERROR: Incorrect format Sigind"
          stop
        endif
        do jj=1, ndim(ii)
          read(fin, *) sigind(1:ndim(ii), jj, ii)
        enddo
        !
        read(fin, *) line
        if (trim(line(1:9))/='Transform') then
          write(*, *) "!!! FATAL ERROR: Incorrect format Transform"
          stop
        endif
        do jj=1, ndim(ii)
          read(fin, *) aa(1:2*ndim(ii))
          do kk=1, ndim(ii)
            Utrans(kk, jj, ii)=aa(2*kk-1)+cmplx_i*aa(2*kk)
          enddo
        enddo
        !
        ! WE ARE IGNORING ALL POSSIBLE "reduced" SUBSPACE CALCULATIONS HERE!!!
        !
        if (mod(ndim(ii)-2, 4).eq.0) then
          corr_l=(ndim(ii)-2)/4
        else
          corr_l=(ndim(ii)-1)/2
        endif
        write(stdout, '(A,1I3)') "     L=", corr_l
        !
        allocate(umat(ndim(ii), ndim(ii)))
        allocate(mtmp(ndim(ii), ndim(ii)))
        allocate(locrot_j(ndim(ii), ndim(ii)))
        ! WIEN2k has its own "local-axis"
        !   which has different definition from us
        !
        call init_symm(symm, locrot, (/0.d0, 0.d0, 0.d0/))
        call inverse_symm(symm)
        locrot_j=cmplx_0
        umat=cmplx_0
        mtmp=cmplx_0
        !
        if (mod(ndim(ii)-2, 4).eq.0) then
          ! ndim(ii)=corr_l*4+2 (SOC case)
          call rotate_J(locrot_j, corr_l, symm)
          call generate_Ylm2C(umat(1:2*corr_l+1, 1:2*corr_l+1), corr_l)
          mtmp(1:2*corr_l+1, 1:2*corr_l+1)=umat(1:2*corr_l+1, 1:2*corr_l+1)
          mtmp(2*corr_l+2:4*corr_l+2, 2*corr_l+2:4*corr_l+2)=umat(1:2*corr_l+1, 1:2*corr_l+1)
        else
          ! ndim(ii)=corr_l*2+1 (non SOC case)
          call rotate_cubic(locrot_j, corr_l, symm)
          call generate_Ylm2C(mtmp, corr_l)
        endif
        !
        call zgemm('C', 'C', ndim(ii), ndim(ii), ndim(ii), &
        cmplx_1, locrot_j, ndim(ii), &
        Utrans(1:ndim(ii), 1:ndim(ii), ii), ndim(ii), &
        cmplx_0, umat, ndim(ii))
        !
        call zgemm('N', 'N', ndim(ii), ndim(ii), ndim(ii), &
        cmplx_1, umat, ndim(ii), &
        mtmp, ndim(ii), &
        cmplx_0, Utrans(1:ndim(ii), 1:ndim(ii), ii), ndim(ii))
        !
        deallocate(umat, mtmp, locrot_j)
        !
      enddo ! imp
      !
      deallocate(aa, locrot)
      !
      close(unit=fin)
      write(stdout, *) " Impurity definitions read."
      !
    endif
    !
    call para_sync_int(basis_map, fulldim)
    call para_sync_int(sigind, maxdim*maxdim*nimp)
    call para_sync_cmplx(Utrans, maxdim*maxdim*nimp)
    !
    if (inode.eq.0) then
      write(stdout, *) " Mapping between impurity problem to full lattice:"
      do ii=1, fulldim
        write(stdout, '(15X,1I5,A,1I5)') ii, '<=>', basis_map(ii)
      enddo
      !
      write(stdout, *) " Transformation Matrices are: "
      do ii=1, nimp
        write(stdout, '(A,1I5)') "    Impurity ", ii
        do jj=1, ndim(ii)
          do kk=1, ndim(ii)
            write(stdout, '(2F9.4,2X)', advance='no') Utrans(jj, kk, ii)
          enddo
          write(stdout, *)
        enddo
      enddo
    endif
    !
    call read_siginp(sigfile, beta)
    !
  end subroutine
  !
  subroutine read_siginp(fname, beta)
    !
    use para,    only: inode, para_sync_int, para_sync_real, para_sync_cmplx
    implicit none
    !
    character(len=80) :: fname
    real(dp)          :: beta
    !
    real(dp), dimension(:), allocatable :: aa
    real(dp) :: beta_sig
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
        omega(ii)=aa(1)
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
    call para_sync_real(omega, nfreq)
    call para_sync_real(sinf, ncol)
    call para_sync_cmplx(sigma, ncol*nfreq)
    !
    if (ismatsubara) then
      beta_sig=(2*nfreq-1.d0)*twopi/(2.d0*omega(nfreq))
      if (abs(beta_sig-beta)>eps4) then
        if (inode.eq.0) then
          write(*, *) "WARNING: beta from input and self energy doesn't match!!!"
        endif
      endif
      !
      deallocate(omega)
    endif
    !
  end subroutine
  !
  subroutine interpolate_sigma(sigpack, w)
    !
    implicit none
    !
    complex(dp), dimension(ncol) :: sigpack
    real(dp) :: w
    !
    integer low, high, mid
    real(dp) :: w1, w2
    !
    if (w>omega(nfreq).or.w<omega(1)) then
      sigpack=cmplx_0
    else
      low=1
      high=nfreq
      do while((high-low)>1)
        mid=(high+low)/2
        if (w>omega(mid)) then
          low=mid
        elseif (w<omega(mid)) then
          high=mid
        else
          high=mid
          low=mid
        endif
      enddo
      !
      if (high.eq.low) then
        sigpack(:)=sigma(:, high)
      else
        w1=omega(low)
        w2=omega(high)
        sigpack(:)=((w-w1)*sigma(:, low)+(w2-w)*sigma(:, high))/(w2-w1)
      endif
      !
    endif
    !
  end subroutine

  subroutine restore_lattice(siglat, sigpack)
    !
    implicit NONE
    !
    complex(dp), dimension(fulldim, fulldim) :: siglat
    complex(dp), dimension(ncol) :: sigpack
    !
    complex(dp), dimension(:, :), allocatable :: sigtmp
    complex(dp), dimension(:, :), allocatable :: mtmp
    !complex(dp), dimension(fulldim, fulldim) :: sigtmp
    !
    integer ii, io1, io2
    siglat=cmplx_0
    !
    do ii=1, nimp
      allocate(sigtmp(ndim(ii), ndim(ii)))
      allocate(mtmp(ndim(ii), ndim(ii)))
      !
      do io1=1, ndim(ii)
        do io2=1, ndim(ii)
          if (sigind(io1, io2, ii)>0) then
            sigtmp(io1, io2)=sigpack(sigind(io1, io2, ii))
          else
            sigtmp(io1, io2)=cmplx_0
          endif
        enddo
      enddo
      !
      call zgemm('N', 'N', ndim(ii), ndim(ii), ndim(ii), &
      cmplx_1, sigtmp, ndim(ii), &
      Utrans(1:ndim(ii), 1:ndim(ii), ii), ndim(ii), &
      cmplx_0, mtmp, ndim(ii))
      !
      call zgemm('C', 'N', ndim(ii), ndim(ii), ndim(ii), &
      cmplx_1, Utrans(1:ndim(ii), 1:ndim(ii), ii), ndim(ii), &
      mtmp, ndim(ii), &
      cmplx_0, sigtmp, ndim(ii))
      !
      do io1=1, ndim(ii)
        do io2=1, ndim(ii)
          siglat(nnlow(ii)+io1-1, nnlow(ii)+io2-1)=sigtmp(io1, io2)
        enddo
      enddo
      !
      deallocate(sigtmp, mtmp)
    enddo
    !
  end subroutine

  subroutine finalize_impurity()
    !
    if (allocated(nnlow)) deallocate(nnlow)
    if (allocated(sigind)) deallocate(sigind)
    if (allocated(Utrans)) deallocate(Utrans)
    if (allocated(basis_map)) deallocate(basis_map)
    if (allocated(sinf)) deallocate(sinf)
    if (allocated(omega)) deallocate(omega)
    if (allocated(sigma)) deallocate(sigma)
    !
  end subroutine

  
END MODULE
