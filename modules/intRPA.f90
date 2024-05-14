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
    integer, dimension(ham%norb) :: mapping
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
    mapping=0
    !
    call para_sync_int0(nFFidx)
    call para_sync_int0(nCCidx)
    !
    if (nCCidx<0) then
      !
      write(*, *) " !!! Incorrect RPA division!"
      stop
      !
    endif
    !
    allocate(FFidx(2, nFFidx))
    if (nCCidx>0) allocate(CCidx(nCCidx))
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
        mapping(blkidx(:))=ii
        !
        do j1=1, blkdim(ii)
          do j2=1, blkdim(ii)
            FFidx(1, jj)=blkidx(j1)
            FFidx(2, jj)=blkidx(j2)
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
        if (mapping(ii).eq.0) then
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
    if (nCCidx>0) call para_sync_int(CCidx, nCCidx)
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
  SUBROUTINE calc_chiRPA(chiff, chicc, chifc, chicf, chi0ff, chi0cc, chi0fc, chi0cf, ff_only, nw)
    !
    use constants,      only : dp, cmplx_0, cmplx_1
    use linalgwrap,     only : sparsemulmat, matmulsparse, invmat
    !
    implicit none
    !
    logical ff_only
    integer nw
    !
    complex(dp), dimension(nFFidx, nFFidx, nw) :: chiff, chi0ff
    complex(dp), dimension(nCCidx, nCCidx, nw) :: chicc, chi0cc
    complex(dp), dimension(nFFidx, nCCidx, nw) :: chifc, chi0fc
    complex(dp), dimension(nCCidx, nFFidx, nw) :: chicf, chi0cf
    !
    complex(dp), dimension(nFFidx, nFFidx) :: Dff, Vff
    complex(dp), dimension(nCCidx, nFFidx) :: tmpcf
    !
    integer ii, iw
    !
    do iw=1, nw
      !
      Dff=cmplx_0
      !
      do ii=1, nFFidx
        Dff(ii, ii)=cmplx_1
      enddo
      !
      ! D_{FF}=(1-\chi^0_{FF}*U_{FF})^{-1}
      ! Vff=Uff*Dff
      !
      call matmulsparse(Dff, chi0ff(:, :, iw),  Uint_cp, idxUcp, nFFidx, nUcp, -1.d0, 1.d0)
      call invmat(Dff, nFFidx)
      call sparsemulmat(Vff, Uint_cp,           Dff,     idxUcp, nFFidx, nUcp,  1.d0, 0.d0)
      !
      ! chiff=Dff*chi0ff
      !
      call zgemm('N', 'N', nFFidx, nFFidx, nFFidx, cmplx_1, Dff, nFFidx, chi0ff(:, :, iw), nFFidx, cmplx_0, chiff(:, :, iw), nFFidx)
      !
      if ((.not. ff_only).and.(nCCidx>0)) then
        !
        ! tmpcf=chi0cf*Vff
        ! chifc=Dff*chi0fc
        ! chicf=chi0cf+tmpcf*chi0ff
        ! chicc=chi0cc+tmpcf*chi0fc
        !
        call zgemm('N', 'N', nCCidx, nFFidx, nFFidx, cmplx_1, chi0cf(:, :, iw), nCCidx, Vff,              nFFidx, cmplx_0, tmpcf,           nCCidx)
        call zgemm('N', 'N', nFFidx, nCCidx, nFFidx, cmplx_1, Dff,              nFFidx, chi0fc(:, :, iw), nFFidx, cmplx_0, chifc(:, :, iw), nFFidx)
        call zgemm('N', 'N', nCCidx, nFFidx, nFFidx, cmplx_1, tmpcf,            nCCidx, chi0ff(:, :, iw), nFFidx, cmplx_1, chicf(:, :, iw), nCCidx)
        call zgemm('N', 'N', nCCidx, nCCidx, nFFidx, cmplx_1, tmpcf,            nCCidx, chi0fc(:, :, iw), nFFidx, cmplx_1, chicc(:, :, iw), nCCidx)
        !
      endif
      !
    enddo ! iw
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
