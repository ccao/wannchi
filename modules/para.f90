MODULE para
  !
  ! This module encapsules the operations with MPI
  !   BRIEFS:
  !     Initialize 
  !         by calling   init_para
  !     Finalize 
  !         by calling   finalize_para
  !     Distribute Calculation
  !         by calling   distribute_calc
  !     Synchronize calculation process
  !         by calling   para_barrier
  !     Synchronize data
  !         by calling   para_sync_TYPE
  !     Summing up data from all process
  !         by calling   para_merge_TYPE
  !     Distribute data to all process
  !         by calling   para_distribute_TYPE
  !     Collect data from all process
  !         by calling   para_collect_TYPE
  !
  use mpi
  !
  implicit none
  !
  integer inode, nnode
  integer first_idx, last_idx
  integer, dimension(:, :), allocatable :: map
  !
  interface para_sync0
    module procedure para_sync_int0, para_sync_real0
  end interface
  !    
  interface para_merge0
    module procedure para_merge_cmplx0, para_merge_real0
  end interface
  !
CONTAINS
  !
SUBROUTINE distribute_calc(nidx)
  !
  implicit none
  !
  integer nidx
  !
#if defined __MPI
  !
  map(:, :)=0
  first_idx=inode*nidx/nnode+1
  last_idx=(inode+1)*nidx/nnode
  map(inode+1, 1)=first_idx-1
  map(inode+1, 2)=last_idx-first_idx+1
  !
  call para_merge_int(map, 2*nnode)
  !
#else
  first_idx=1
  last_idx=nidx
#endif
  !
END SUBROUTINE

SUBROUTINE init_para(codename)
  !
  use constants
  !
  implicit none
  !
  character(*) codename
  !
#if defined __MPI
  integer ierr
  !
  CALL mpi_init(ierr)
  !
  CALL mpi_comm_rank(mpi_comm_world, inode, ierr)
  CALL mpi_comm_size(mpi_comm_world, nnode, ierr)
  !
  if (inode.eq.0) write(stdout, *) trim(codename)//" running on ", nnode, " nodes..."
  !
  allocate(map(nnode, 2))
  !
#else
  inode=0
  nnode=1
  write(stdout, *) trim(codename)//" serial ..."
#endif
  !
END SUBROUTINE

SUBROUTINE para_barrier()
  !
  implicit none
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_barrier(mpi_comm_world, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_int0(dat)
  !
  implicit none
  !
  integer :: dat
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_bcast(dat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_int(dat, dat_size)
  !
  implicit none
  !
  integer :: dat(*)
  integer :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_bcast(dat, dat_size, MPI_INTEGER, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real(dat, dat_size)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) :: dat(*)
  integer  :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real0(dat)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) :: dat
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_bcast(dat, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_cmplx(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  complex(dp) :: dat(*)
  integer :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_int(dat, dat_size)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer  :: dat(*)
  integer  :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_allreduce(MPI_IN_PLACE, dat, dat_size, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real0(dat)
  !
  use constants, only : dp
  !
  implicit none
  !
  real(dp) :: dat
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_allreduce(MPI_IN_PLACE, dat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real(dat, dat_size)
  !
  use constants, only :dp
  !
  implicit none
  !
  real(dp) :: dat(*)
  integer  :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_allreduce(MPI_IN_PLACE, dat, dat_size, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx0(dat)
  !
  use constants, only : dp
  !
  implicit none
  !
  complex(dp) :: dat
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_allreduce(MPI_IN_PLACE, dat, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  complex(dp) :: dat(*)
  integer :: dat_size
  !
#if defined __MPI
  !
  integer ierr
  CALL mpi_allreduce(MPI_IN_PLACE, dat, dat_size, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_collect_int(fulldat, dat)
  !
  use constants, only : dp
  implicit none
  !
  integer :: fulldat(*)
  complex(dp) :: dat(*)
  !
#if defined __MPI
  !
  ! We need to use mpi_type
  integer ierr
  !
  call mpi_gatherv(dat, map(inode+1, 2), MPI_INTEGER, fulldat, map(:, 2), map(:, 1), MPI_INTEGER, 0, mpi_comm_world, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE para_collect_cmplx(fulldat, dat, blk_size)
  !
  use constants, only : dp
  implicit none
  !
  complex(dp) :: fulldat(*)
  complex(dp) :: dat(*)
  integer :: blk_size
  !
#if defined __MPI
  !
  ! We need to use mpi_type
  integer ierr, blk_cmplx
  !
  call mpi_type_contiguous(blk_size, MPI_DOUBLE_COMPLEX, blk_cmplx, ierr)
  ! new type object created
  call mpi_type_commit(blk_cmplx, ierr)
  ! now type1 can be used for communication
  call mpi_gatherv(dat, map(inode+1, 2), blk_cmplx, fulldat, map(:, 2), map(:, 1), blk_cmplx, 0, mpi_comm_world, ierr)
  call mpi_type_free(blk_cmplx, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE para_distribute_cmplx(fulldat, dat, blk_size)
  !
  use constants, only : dp
  implicit none
  !
  complex(dp) :: fulldat(*)
  complex(dp) :: dat(*)
  integer :: blk_size
  !
#if defined __MPI
  !
  ! We need to use mpi_type
  integer ierr, blk_cmplx
  !
  call mpi_type_contiguous(blk_size, MPI_DOUBLE_COMPLEX, blk_cmplx, ierr)
  ! new type object created
  call mpi_type_commit(blk_cmplx, ierr)
  ! now type1 can be used for communication
  CALL mpi_scatterv(fulldat, map(:, 2), map(:, 1), blk_cmplx, dat, map(inode+1, 2), blk_cmplx, 0, mpi_comm_world, ierr)
  call mpi_type_free(blk_cmplx, ierr)
  !
#endif
  !
END SUBROUTINE

SUBROUTINE finalize_para
  !
  implicit none
  !
#if defined __MPI
  integer ierr
  CALL mpi_finalize(ierr)
#endif
  !
END SUBROUTINE

END MODULE
