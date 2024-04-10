MODULE lattice
  !
  ! This module encapsules the operations of lattice problem
  !   BRIEFS:
  !     Initialize lattice structure information
  !         by calling   read_posfile
  !
  use constants,   only : dp
  use wanndata,    only : wannham
  use simp_module, only : simp
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
  integer nimp
  ! Number of impurities
  !
  TYPE(simp), dimension(:), allocatable  :: imp
  ! Impurity problems 
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
  complex(dp), dimension(:), allocatable :: omega
  ! Frequency mesh
  !
  complex(dp), dimension(:), allocatable :: sinf
  ! Sigma(inf)
  !
  complex(dp), dimension(:, :), allocatable :: sigpack
  ! Sigma(w) in packed form
  !
  integer, dimension(:), allocatable     :: partition
  ! Partition of basis
  !  partition(ham%norb), partition(c)=0, partition(f_i)=i
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
    use simp_module,  only : finalize_simp
    !
    implicit none
    !
    integer ii
    !
    call finalize_wann(ham, .true.)
    !
    do ii=1, nimp
      call finalize_simp(imp(ii))
    enddo
    !
    if (allocated(partition)) deallocate(partition)
    if (allocated(sinf)) deallocate(sinf)
    if (allocated(sigpack)) deallocate(sigpack)
    !
  end subroutine
  !
  subroutine read_sigfile(fn)
    !
    ! Soo ...
    ! Vdc ...
    ! S(w_1)
    ! ...
    !
    use constants,    only : stdout, dp, fin, cmplx_i
    use para,         only : inode, para_sync_int, para_sync_cmplx
    !
    implicit none
    !
    character(*) fn
    !
    integer, dimension(2) :: tt
    integer ii, jj
    !
    real(dp), dimension(:), allocatable :: rr
    !
    if (inode.eq.0) then 
      open(unit=fin, file=trim(fn))
      read(fin, *) tt        ! nfreq  nbath
    endif
    !
    call para_sync_int(tt, 2)
    nw=tt(1)
    nbath=tt(2)
    !
    allocate(omega(nw))
    allocate(sinf(nbath), sigpack(nbath, nw))
    allocate(rr(2*nbath+2))
    !
    if (inode.eq.0) then
      !
      write(stdout, '(A,1I3,A,1I6,A)') "   # Reading self-energy of ", nbath, " baths with ", nw, " frequencies."
      !
      read(fin, *) sinf                 ! Soo
      read(fin, *) rr(1:nbath)          ! Vdc
      do ii=1, nbath
        sinf(ii)=sinf(ii)-rr(ii)
      enddo
      !
      write(stdout, *) " # Sinf : "
      do ii=1, nbath
        write(stdout, '(2F14.9,2X)', advance='no') sinf(ii)
      enddo
      write(stdout, *)
      !
      do jj=1, nw
        read(fin, *) rr      ! s(jj)
        omega(jj)=rr(1)+rr(2)*cmplx_i
        do ii=1, nbath
          sigpack(ii, jj)=rr(2*ii+1)+rr(2*ii+2)*cmplx_i
        enddo
      enddo
      !
      close(fin)
      !
    endif
    !
    deallocate(rr)
    !
    call para_sync_cmplx(omega, nw)
    call para_sync_cmplx(sinf, nbath)
    call para_sync_cmplx(sigpack, nbath*nw)
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
  subroutine read_impfile(fn)
    !
    !  Line 1 : nimp
    !  Line 2 :  ndim1, l1
    !  Line 3 :  gidx1 gidx2 gidx3 ...
    !  Line 4.. :  sigidx11 sigidx12
    !
    use constants,    only : dp, fin, cmplx_i, stdout
    use para,         only : inode, para_sync0, para_sync_int, para_sync_cmplx
    use simp_module,  only : init_simp, simp
    !
    implicit none
    !
    character(*) fn
    !
    integer, dimension(2) :: tt
    integer ii, jj, kk
    real(dp), dimension(:), allocatable :: aa
    !
    if (inode.eq.0) then
      open(unit=fin, file=trim(fn))
      read(fin, *) nimp
      !
      write(stdout, '(A,1I5,A)') "  # Reading definition of ", nimp, " impurities."
      !
    endif
    !
    call para_sync0(nimp)
    allocate(imp(nimp))
    !
    do ii=1, nimp
      !
      if (inode.eq.0) then
        read(fin, *)           !  IMP idx
        read(fin, *) tt        !  ndim,  l
        write(stdout, '(A,1I3,A,1I2,A,1I2,A,1I1)') "     # impurity ", ii, ": ", tt(1), "x", tt(1), " l= ", tt(2)  
      endif
      !
      call para_sync_int(tt, 2)
      !
      call init_simp(imp(ii), tt(1), tt(2))
      !
      if (inode.eq.0) then
        read(fin, *)           ! BasisMap
        read(fin, *) imp(ii)%gidx(:)
        !
        write(stdout, *) "   #  Mapping from impurity to global: "
        do jj=1, imp(ii)%ndim
          write(stdout, '(1I5,A,1I5)') jj, " <==> ", imp(ii)%gidx(jj)
        enddo
        !
        read(fin, *)           ! Sigidx
        do jj=1, imp(ii)%ndim
          read(fin, *) imp(ii)%sigidx(jj, :)
        enddo
        !
        write(stdout, *) "   #  Self energy structure: "
        do jj=1, imp(ii)%ndim
          do kk=1, imp(ii)%ndim
            write(stdout, '(1I5)', advance='no') imp(ii)%sigidx(jj, kk)
          enddo
          write(stdout, *)
        enddo
        !
        read(fin, *)           ! Transformation Matrix
        !
        allocate(aa(2*imp(ii)%ndim))
        do jj=1, imp(ii)%ndim
          read(fin, *) aa
          do kk=1, imp(ii)%ndim
            imp(ii)%Utrans(kk, jj)=aa(kk*2-1)+aa(kk*2)*cmplx_i
          enddo
        enddo
        !
        deallocate(aa)
        !
      endif
      !
      call para_sync_int(imp(ii)%gidx, imp(ii)%ndim)
      call para_sync_int(imp(ii)%sigidx, imp(ii)%ndim*imp(ii)%ndim)
      call para_sync_cmplx(imp(ii)%Utrans, imp(ii)%ndim*imp(ii)%ndim)
      !
    enddo
    !
    if (inode.eq.0) then
      close(unit=fin)
      write(stdout, *) " #  Done."
    endif
    !
  end subroutine
  !
  subroutine lattice_reduce(sig, sigmat, norb)
    !
    use constants,     only : dp, cmplx_0
    use simp_module,   only : simp, matrix2pack
    !
    integer                                    :: norb
    complex(dp), dimension(norb, norb)         :: sigmat
    complex(dp), dimension(nbath)              :: sig
    !
    integer, dimension(nbath)                  :: ndegen
    complex(dp), dimension(:, :), allocatable  :: sigimp
    !
    integer ii, io, jo
    !
    sig=cmplx_0
    ndegen=0
    !
    do ii=1, nimp
      !
      allocate(sigimp(imp(ii)%ndim, imp(ii)%ndim))
      !
      do io=1, imp(ii)%ndim
        do jo=1, imp(ii)%ndim
          sigimp(io, jo)=sigmat(imp(ii)%gidx(io), imp(ii)%gidx(jo))
          if (imp(ii)%sigidx(io, jo)>0) then
            ndegen(imp(ii)%sigidx(io, jo))=ndegen(imp(ii)%sigidx(io, jo))+1
          endif
        enddo
      enddo
      !
      CALL matrix2pack(sig, sigimp, imp(ii))
      !
      deallocate(sigimp)
      !
    enddo
    !
    do ii=1, nbath
      sig(ii)=sig(ii)/ndegen(ii)
    enddo
    !
  end subroutine
  !
  subroutine restore_lattice(sigmat, norb, sig)
    !
    use constants,     only : dp, cmplx_0, stdout
    use simp_module,   only : simp, restore_matrix
    use para,          only : inode
    !
    implicit none
    !
    integer                                    :: norb
    complex(dp), dimension(norb, norb)         :: sigmat
    complex(dp), dimension(nbath)              :: sig
    !
    complex(dp), dimension(:, :), allocatable  :: sigimp
    !
    integer ii, io, jo
    !
    sigmat=cmplx_0
    !
    do ii=1, nimp
      !
      allocate(sigimp(imp(ii)%ndim, imp(ii)%ndim))
      !
      CALL restore_matrix(sigimp, sig, imp(ii))
      !
      do io=1, imp(ii)%ndim
        do jo=1, imp(ii)%ndim
          sigmat(imp(ii)%gidx(io), imp(ii)%gidx(jo))=sigimp(io, jo)
        enddo
      enddo
      !
      deallocate(sigimp)
      !
    enddo
    !
  end subroutine
  !
  subroutine print_impurity(mat)
    !
    use constants,     only : dp, stdout
    use para,          only : inode
    !
    implicit none
    !
    complex(dp), dimension(:, :) :: mat
    integer                      :: ii, i1, i2
    !
    do ii=1, nimp
      !
      write(stdout, '(A,1I3)') "    # Impurity ", ii
      !
      do i1=1, imp(ii)%ndim
        do i2=1, imp(ii)%ndim
          write(stdout, '(2F9.4,2X)', advance='no') mat(imp(ii)%gidx(i1), imp(ii)%gidx(i2))
        enddo
        write(stdout, *)
      enddo
      !
    enddo
    !
  end subroutine
  !
  subroutine interpolate_single_sigma(sig, z)
    !
    use constants, only : dp, twopi, eps6, cmplx_i, eps12
    !
    implicit none
    !
    complex(dp), dimension(nbath) :: sig
    complex(dp) :: z
    !
    complex(dp), dimension(nbath) :: sigtmp
    real(dp) :: w  ! The pure imaginary part of z
    !
    integer low, high, mid
    real(dp) :: w1, w2, s11, s12, s21, s22
    !
    real(dp) :: ff
    real(dp) :: a, b, c, d
    integer ii
    !
    w=aimag(z)
    beta=(2*nw-1.d0)*twopi/(2.d0*aimag(omega(nw)))
    ff=(abs(w)*beta/twopi+0.5d0)
    ii=nint(ff)
    if (abs(ff-ii)<eps6 .and. ii<=nw) then
      ! Exactly on the mesh
      !
      sigtmp(:)=sigpack(:, ii)
      !
    else if (ff>nw) then
      ! Out of the mesh
      !
      do ii=1, nbath
        ! Find out parameters generating HB-1 tail
        !
        w1=aimag(omega(nw-10))
        w2=aimag(omega(nw))
        !
        s11=real(sigpack(ii, nw-10))
        s12=aimag(sigpack(ii, nw-10))
        s21=real(sigpack(ii, nw))
        s22=aimag(sigpack(ii, nw))
        !
        a=(s21-s11)/(1.d0/(w2*w2)-1.d0/(w1*w1))
        d=(s21*w2*w2-s11*w1*w1)/(w2*w2-w1*w1)
        !
        if (abs(s22)<eps6 .or. abs(s12)<eps6) then
          b=0.d0
          c=0.d0
        else
          b=-(w2*s22-w1*s12)/(s22/w2-s12/w1)
          c=(w2*w2-w1*w1)/(w2/s22-w1/s12)
        endif
        ! Calculate the correct tail value
        !
        sigtmp(ii)=d+a/(w*w)+cmplx_i*c*abs(w)/(w*w+b)
        !
      enddo
      !
    else if (ff>ii .and. ii<nw) then
      ! Between ii and ii+1
      !
      sigtmp(:)=(ii+1-ff)*sigpack(:, ii)+(ff-ii)*sigpack(:, ii+1)
      !
    else if (ff<ii .and. ii>0) then
      ! Between ii and ii-1
      !
      sigtmp(:)=(ii-ff)*sigpack(:, ii-1)+(ff-ii+1)*sigpack(:, ii)
      !
    else
      write(*, *)  '!!! FATAL: Incorrect matsubara frequency! Are you sure?'
      stop
    endif
    !
    if (w>0) then
      sig(:)=sigtmp(:)
    else
      sig(:)=conjg(sigtmp(:))
    endif
    !
  end subroutine
  !
  subroutine interpolate_sigmas(sig_int, w_int, n_int)
    !
    use constants, only : dp
    !
    implicit none
    !
    integer        :: n_int
    complex(dp), dimension(nbath, n_int) :: sig_int
    complex(dp), dimension(n_int)        :: w_int
    !
    integer ii
    do ii=1, n_int
      CALL interpolate_single_sigma(sig_int(:, ii), w_int(ii))
    enddo
    !
  end subroutine
  !
END MODULE
