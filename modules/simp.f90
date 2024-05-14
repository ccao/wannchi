!
!   simp.f90
!   
!
!   Created by Chao Cao on 06/01/23.
!   Copyright 2023 CC. All rights reserved.
!
MODULE simp_module
  !
  ! This module encapsules the operations on Wannier Hamiltonian
  !   BRIEFS:
  !     Initialize 
  !         by calling   init_simp
  !     Finalize
  !         by calling   finalize_simp
  !     Transform from compact form to (local) matrix form
  !         by calling   restore_matrix
  !
  use constants,   only : dp
  !
  implicit none
  !
  TYPE simp
    !
    integer ndim
    ! Dimension of this impurity
    !
    integer lang
    ! angular momentum number of this impurity
    !
    integer, dimension(:), allocatable :: gidx
    ! Global orbital indices : gidx(ndim)
    !
    complex(dp), dimension(:, :), allocatable :: Utrans
    ! Transform matrix from DMFT basis to cubic harmonix
    !
    integer, dimension(:, :), allocatable :: sigidx
    ! Self energy mapping: Sigma(i, j)=sigpack(sigidx(i,j))
    !
  END TYPE
  !
  contains
  !
  subroutine init_simp(imp, ndim, l)
    !
    implicit none
    !
    TYPE(simp), intent(out)  :: imp
    integer, intent(in)      :: ndim, l
    !
    imp%ndim=ndim
    imp%lang=l
    !
    allocate(imp%gidx(ndim))
    allocate(imp%Utrans(ndim, ndim))
    allocate(imp%sigidx(ndim, ndim))
    !
  end subroutine
  !
  subroutine finalize_simp(imp)
    !
    implicit none
    !
    TYPE(simp), intent(inout)  :: imp
    !
    if (allocated(imp%gidx))   deallocate(imp%gidx)
    if (allocated(imp%Utrans)) deallocate(imp%Utrans)
    if (allocated(imp%sigidx)) deallocate(imp%sigidx)
    !
  end subroutine
  !
  subroutine matrix2pack(sigpack, sigmat, imp)
    !
    use constants,     only : dp, cmplx_0, cmplx_1
    !
    implicit none
    !
    TYPE(simp), intent(in)                                  :: imp
    complex(dp), dimension(imp%ndim, imp%ndim), intent(in)  :: sigmat
    complex(dp), dimension(:), intent(out)                  :: sigpack
    !
    complex(dp), dimension(imp%ndim, imp%ndim) :: stmp1, stmp2
    integer ii, jj
    !
    call zgemm('N', 'C', imp%ndim, imp%ndim, imp%ndim, &
                cmplx_1, sigmat, imp%ndim, &
                imp%Utrans, imp%ndim, &
                cmplx_0, stmp1, imp%ndim)
    !
    call zgemm('N', 'N', imp%ndim, imp%ndim, imp%ndim, &
                cmplx_1, imp%Utrans, imp%ndim, &
                stmp1, imp%ndim, &
                cmplx_0, stmp2, imp%ndim)
    !
    do ii=1, imp%ndim
      do jj=1, imp%ndim
        !
        if (imp%sigidx(ii, jj)>0) then
          sigpack(imp%sigidx(ii, jj))=sigpack(imp%sigidx(ii, jj))+sigmat(ii, jj)
        endif
        !
      enddo
    enddo
    !
  end subroutine
  !
  subroutine restore_matrix(sigmat, sigpack, imp)
    !
    use constants,      only : dp, cmplx_0, cmplx_1
    !
    implicit none
    !
    TYPE(simp), intent(in)                                  :: imp
    complex(dp), dimension(imp%ndim, imp%ndim), intent(out) :: sigmat
    complex(dp), dimension(:), intent(in)                   :: sigpack
    !
    complex(dp), dimension(imp%ndim, imp%ndim) :: stmp
    integer ii, jj
    !
    sigmat=cmplx_0
    !
    do ii=1, imp%ndim
      do jj=1, imp%ndim
        !
        if (imp%sigidx(ii, jj)>0) then
          sigmat(ii, jj)=sigpack(imp%sigidx(ii, jj))
        endif
        !
      enddo
    enddo
    !
    call zgemm('N', 'N', imp%ndim, imp%ndim, imp%ndim, &
                cmplx_1, sigmat, imp%ndim, &
                imp%Utrans, imp%ndim, &
                cmplx_0, stmp, imp%ndim)
    !
    call zgemm('C', 'N', imp%ndim, imp%ndim, imp%ndim, &
                cmplx_1, imp%Utrans, imp%ndim, &
                stmp, imp%ndim, &
                cmplx_0, sigmat, imp%ndim)
    !
  end subroutine
  !
END MODULE