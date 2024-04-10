MODULE IntRPA
  !
  use constants,   only :dp
  !
  implicit none
  !
  integer nFFidx
  ! matrix dimension for FF-block
  !
  integer nCCidx
  ! matrix dimension for CC-block
  !
  integer, dimension(:, :), allocatable :: FFidx
  ! indice pair mapping between F and ff in full matrix
  !   FFidx(1, i)=i1, FFidx(2, i)=i2
  !
  integer, dimension(:), allocatable :: CCidx
  ! indice mapping in CC block
  !   CCidx(i)=
  !
  integer nUcp
  ! Number of U matrix elements in compact form
  real(dp), dimension(:), allocatable :: Uint_cp
  ! U matrix in compact form
  integer, dimension(:, :), allocatable :: idxUcp
  ! nonzero element indices in full matrix for U matrix
  ! U_{FF}(idxUcp(1,ii), idxUcp(2,ii))=UFF_{cp}(ii)
  !
 contains
  !
  SUBROUTINE read_RPA
    !
    use constants,    only : dp, fin
    use lattice,      only : ham
    use para,         only : inode, para_sync_int0, para_sync_int, para_sync_real
    !
    implicit none
    !
    integer nffblk
    integer, dimension(:), allocatable :: blkdim
    integer, dimension(:), allocatable :: blkidx
    integer, dimension(:), allocatable :: partition
    integer ii, jj, j1, j2, i1, i2
    !
    ! Figure out dimensions...
    !
    if (inode.eq.0) then
      !
      open(unit=fin, file="RPA.inp")
      read(fin, *) nffblk
      allocate(blkdim(nffblk))
      read(fin, *) blkdim
      !
      nFFidx=0
      nCCidx=ham%norb
      !
      do ii=1, nffblk
        nFFidx=nFFidx+blkdim(ii)*blkdim(ii)
        nCCidx=nCCidx-blkdim(ii)
      enddo
      !
    endif
    !
    call para_sync_int0(nFFidx)
    call para_sync_int0(nCCidx)
    !
    allocate(FFidx(2, nFFidx), CCidx(nCCidx))
    allocate(partition(ham%norb))
    partition(:)=0
    !
    ! Construct FFblock / CCblock
    !
    if (inode.eq.0) then
      !
      jj=1
      do ii=1, nffblk
        !
        allocate(blkidx(blkdim(ii)))
        !
        read(fin, *) blkidx(:)
        partition(blkidx(:))=ii
        !
        do j1=1, blkdim(ii)
          do j2=1, blkdim(ii)
            FFidx(1, jj)=j1
            FFidx(2, jj)=j2
            jj=jj+1
          enddo
        enddo
        !
        deallocate(blkidx)
        !
      enddo
      !
      jj=1
      do ii=1, ham%norb
        !
        if (partition(ii).eq.0) then
          CCidx(jj)=ii
          jj=jj+1
        endif
        !
      enddo
      !
      read(fin, *) nUcp
      !
    endif
    !
    call para_sync_int(FFidx, nFFidx*2)
    call para_sync_int(CCidx, nCCidx)
    call para_sync_int0(nUcp)
    !
    if (nUcp.ne.0) then
      !
      allocate(Uint_cp(nUcp), idxUcp(2, nUcp))
      !
      if (inode.eq.0) then
        !
        do ii=1, nUcp
          !
          read(fin, *) i1, i2, j1, j2, Uint_cp(ii)
          call find_ffidx(jj, i1, i2)
          idxUcp(1, ii)=jj
          call find_ffidx(jj, j1, j2)
          idxUcp(2, ii)=jj
          !
        enddo
        !
      endif
      !
      call para_sync_int(idxUcp, nUcp*2)
      call para_sync_real(Uint_cp, nUcp)
      !
    endif
    !
    if (inode.eq.0) then
      !
      deallocate(blkdim)
      close(fin)
      !
    endif
    !
  END SUBROUTINE
  !
  SUBROUTINE find_ffidx(ii, i1, i2)
    !
    implicit none
    !
    integer ii, i1, i2
    integer jj
    !
    ii=0
    !
    do jj=1, nFFidx
      !
      if (FFidx(1, jj)==i1 .and. FFidx(2, jj)==i2) then
        ii=jj
        exit
      endif
      !
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_RPA
    !
    implicit none
    !
    if (allocated(FFidx))   deallocate(FFidx)
    if (allocated(CCidx))   deallocate(CCidx)
    if (allocated(Uint_cp)) deallocate(Uint_cp)
    if (allocated(idxUcp))  deallocate(idxUcp)
    !
  END SUBROUTINE
  !
END MODULE
