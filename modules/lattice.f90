MODULE lattice
  !
  ! This module encapsules the operations of lattice problem
  !   BRIEFS:
  !     Initialize lattice structure information
  !         by calling   read_posfile
  !
  use constants,   only : dp
  use wanndata,    only : wannham
  !
  implicit none
  !
  real(dp), dimension(3, 3) :: avec
  ! Lattice vectors
  !  a1=avec(:, 1); a2=avec(:, 2); a3=avec(:, 3)
  !
  real(dp), dimension(3, 3) :: bvec
  ! Reciprocal lattice vectors
  !  b1=bvec(1, :); b2=bvec(2, :); b3=bvec(3, :) 
  !
  integer nsite
  ! Number of sites in the lattice
  !
  integer, dimension(:), allocatable     :: zat
  ! Atomic Numbers
  !
  real(dp), dimension(:, :), allocatable :: xat
  ! Fractional coordinates of atoms
  !
  integer, dimension(:), allocatable     :: nbasis
  ! Number of basis at each site
  !
  real(dp)  :: nelec
  ! Total number of electrons (valence)
  !
  TYPE(wannham)                          :: ham
  ! Hamiltonian associated with this lattice
  !
  integer nbath
  ! Number of bath functions
  !
  real(dp) beta
  ! Temperature
  !
  integer nw
  ! Number of frequencies in frequency mesh
  !
  integer nk1, nk2, nk3
  ! BZ mesh size
  !
  integer nkirr
  ! Number of irreducible Kpts in BZ
  !
  real(dp), dimension(:), allocatable    :: kwt
  ! Weight of Kpts in BZ
  !
  real(dp), dimension(:, :), allocatable :: kvec
  ! Kvectors
  !
  contains
  !
  subroutine finalize_lattice_structure()
    !
    implicit none
    !
    if (allocated(zat)) deallocate(zat)
    if (allocated(xat)) deallocate(xat)
    if (allocated(nbasis)) deallocate(nbasis)
    !
  end subroutine
  !
  subroutine finalize_lattice_kmesh()
    !
    implicit none
    !
    if (allocated(kwt)) deallocate(kwt)
    if (allocated(kvec)) deallocate(kvec)
    !
  end subroutine
  !
  subroutine finalize_lattice_ham()
    !
    use wanndata,     only : finalize_wann
    !
    implicit none
    !
    integer ii
    !
    call finalize_wann(ham, .true.)
    !
  end subroutine
  !
  subroutine read_posfile(fn)
    !
    ! READ posfile
    !  Line 1 : comment
    !  Line 2 :  universal scaling factor
    !  Line 3 :   A1
    !  Line 4 :   A2
    !  Line 5 :   A3
    !  Line 6 :  nsite soc
    !  Line 7..nsite+6 :  Zat   X1  X2  X3  Nbasis
    !
    use constants,    only : dp, fin
    use para,         only : inode, para_sync_int0, para_sync_int, para_sync_real
    use linalgwrap,   only : invmat
    !
    implicit none
    !
    character(*), intent(in) :: fn
    !
    real(dp)  :: alat
    real(dp), dimension(3)  :: xx
    integer   :: ii, jj, nn, ispin
    !
    if (inode.eq.0) then
      open(unit=fin, file=trim(fn))
      !
      read(fin, *) ! First line is comment
      !
      read(fin, *) alat
      !
      do ii=1, 3
        read(fin, *) xx
        avec(:, ii)=xx(:)*alat
      enddo
      !
      bvec(:, :)= avec(:, :)
      !
      call invmat(bvec, 3)
      !
      read(fin, *) nsite, ispin
      !
    endif
    !
    call para_sync_real(avec, 9)
    call para_sync_real(bvec, 9)
    !
    call para_sync_int0(nsite)
    call para_sync_int0(ispin)
    !
    allocate(xat(3, nsite))
    allocate(zat(nsite))
    allocate(nbasis(nsite))
    !
    if (inode.eq.0) then
      !
      do ii=1, nsite
        read(fin, *) zat(ii), xat(:, ii), nbasis(ii)
      enddo
      !
      close(unit=fin)
      !
    endif
    !
    call para_sync_int(zat, nsite)
    call para_sync_int(nbasis, nsite)
    call para_sync_real(xat, 3*nsite)
    !
    nn=1
    do ii=1, nsite
      do jj=1, nbasis(ii)
        ham%tau(:, nn)=xat(:, ii)
        nn=nn+1
      enddo
    enddo
    !
    if (ispin>0) then
      do ii=1, nsite
        do jj=1, nbasis(ii)
          ham%tau(:, nn)=xat(:, ii)
          nn=nn+1
        enddo
      enddo
    endif
    !
  end subroutine
  !
  subroutine read_kmesh(fn)
    !
    use constants,     only : dp, fin, stdout
    use para,          only : inode, para_sync_int, para_sync_real
    !
    implicit none
    !
    character(*)    :: fn
    !
    integer ik1, ik2, ik3, iik
    integer, dimension(4) :: tt
    !
    if (inode.eq.0) then
      !
      open(unit=fin, file=trim(fn))
      !
      read(fin, *)      ! COMMENT
      read(fin, *) iik  ! SWITCH
      !
      if (iik.eq.0) then
        ! AUTOMATIC K_MESH
        read(fin, *)    ! Always use Gamma-centered
        read(fin, *) nk1, nk2, nk3
        !
        nkirr=nk1*nk2*nk3
      else
        read(fin, *)
        !
        nkirr=iik
      endif
      !
      tt(1)=nk1
      tt(2)=nk2
      tt(3)=nk3 
      tt(4)=nkirr
      !
    endif
    !
    CALL para_sync_int(tt, 4)
    !
    nk1=tt(1)
    nk2=tt(2)
    nk3=tt(3)
    nkirr=tt(4)
    !
    allocate(kwt(nkirr), kvec(3, nkirr))
    !
    if (inode.eq.0) then
      if (iik.eq.0) then
        !
        do ik1=0, nk1-1
          do ik2=0, nk2-1
            do ik3=0, nk3-1
              iik=ik1*nk2*nk3+ik2*nk3+ik3+1
              kwt(iik)=1.d0
              kvec(1, iik)=ik1*1.d0/nk1
              kvec(2, iik)=ik2*1.d0/nk2
              kvec(3, iik)=ik3*1.d0/nk3
            enddo
          enddo
        enddo
        !
        iik=0
        !
      else
        ! IBZKPT form
        do iik=1, nkirr
          read(fin, *) kvec(:, iik), ik1
          kwt(iik)=ik1*1.d0
        enddo
        !
      endif
      !
      close(unit=fin)
      !
    endif
    !
    call para_sync_real(kwt, nkirr)
    call para_sync_real(kvec, nkirr*3)
    !
    kwt(:)=kwt(:)/SUM(kwt(:))
    !
    if (inode.eq.0) then
      !
      if (iik==0) then
        write(stdout, '(A)', advance='no') " Automatically Generated K-mesh of "
        write(stdout, '(1I5,1A1,1I5,1A1,1I5)') nk1, 'x', nk2, 'x', nk3
      else
        write(stdout, '(A,1I5,A)') " Specified K-mesh of ", nkirr, " K-points"
      endif
    endif
    !
  end subroutine
  !
END MODULE
