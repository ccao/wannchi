MODULE para
  !
  !use mpi
  !
  implicit none
  !
  integer inode, nnode
  integer first_k, last_k
  !
  interface para_merge0
    module procedure para_merge_cmplx0
  end interface
  !
!INTERFACE para_sync
!  MODULE PROCEDURE para_sync_int0, para_sync_real0, para_sync_real1, para_sync_cmplx1, para_sync_real2, para_sync_cmplx3
!END INTERFACE
  !
!INTERFACE para_merge
!  MODULE PROCEDURE para_merge_real0, para_merge_real1, para_merge_real2, para_merge_cmplx0, para_merge_cmplx1, &
!    para_merge_cmplx3, para_merge_cmplx4
!END INTERFACE
  !
CONTAINS
  !
SUBROUTINE distribute_k()
  !
  use banddata, only : nkpt
  !
  implicit none
  !
#if defined __MPI
  first_k=inode*nkpt/nnode+1
  last_k=(inode+1)*nkpt/nnode
#else
  first_k=1
  last_k=nkpt
#endif
  !
END SUBROUTINE

SUBROUTINE init_para()
  !
  use constants
  !
  implicit none
  !
#if defined __MPI
  integer ierr
  !
  CALL mpi_init(ierr)
  !
  CALL mpi_comm_rank(mpi_comm_world, inode, ierr)
  CALL mpi_comm_size(mpi_comm_world, nnode, ierr)
  !
  if (inode.eq.0) write(stdout, *) "# WannChi running on ", nnode, " nodes..."
#else
  inode=0
  nnode=1
  write(stdout, *) "# WannChi serial ..."
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
