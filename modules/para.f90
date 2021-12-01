MODULE para
  !
  use mpi
  !
  implicit none
  !
  integer inode, nnode
  integer first_idx, last_idx
  integer, dimension(:, :), allocatable :: map
  !
  interface para_merge0
    module procedure para_merge_cmplx0
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
  integer ierr
  integer tmap(nnode, 2)
  tmap(:, :)=map(:, :)*blk_size
  call mpi_gatherv(dat, tmap(inode+1, 2), MPI_DOUBLE_COMPLEX, fulldat, tmap(:, 2), tmap(:, 1), MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
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
  integer ierr
  integer tmap(nnode, 2)
  tmap(:, :)=map(:, :)*blk_size
  CALL mpi_scatterv(fulldat, tmap(:, 2), tmap(:, 1), MPI_DOUBLE_COMPLEX, dat, tmap(inode+1, 2), MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
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
