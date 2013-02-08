MODULE para
  !
#if defined __MPI
  use mpi
#endif
  !
  implicit none
  !
  integer inode, nnode
  !
INTERFACE para_sync
  MODULE PROCEDURE para_sync_int0, para_sync_real0, para_sync_real1, para_sync_real2, para_sync_cmplx3
END INTERFACE
  !
INTERFACE para_merge
  MODULE PROCEDURE para_merge_real1, para_merge_real2, para_merge_cmplx1, para_merge_cmplx3
END INTERFACE
  !
CONTAINS
  !
SUBROUTINE init_para()
  !
  use constants
  !
  implicit none
  !
  integer ierr
  !
#if defined __MPI
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

SUBROUTINE para_sync_int0(dat)
  !
  implicit none
  !
  integer :: dat
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
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
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real1(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real2(dat, size1, size2)
  !
  use constants, only: dp
  !
  implicit none
  !
  integer :: size1, size2
  real(dp) :: dat(1:size1, 1:size2)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_cmplx3(dat, size1, size2, size3)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: size1, size2, size3
  complex(dp) :: dat(1:size1, 1:size2, 1:size3)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2*size3
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE
  
SUBROUTINE para_merge_real1(dat, dat_size)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
  integer loc_size
  real(dp), allocatable :: buffer(:)
  !
#if defined __MPI
  loc_size=dat_size/nnode
  allocate(buffer(1:loc_size))
  !
  buffer(:)=dat(inode*loc_size+1:(inode+1)*loc_size)
  !
  CALL mpi_allgather(buffer, loc_size, MPI_DOUBLE_PRECISION, dat, loc_size, MPI_DOUBLE_PRECISION, mpi_comm_world, ierr)
  !
  deallocate(buffer)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real2(dat, size1, size2)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer :: size1, size2
  real(dp) :: dat(1:size1, size2)
  !
  integer ierr
  !
  integer loc_size
  real(dp), allocatable :: buffer(:, :)
  !
#if defined __MPI
  loc_size=size2/nnode
  allocate(buffer(1:size1, 1:loc_size))
  !
  buffer(:,:)=dat(:, inode*loc_size+1:(inode+1)*loc_size)
  !
  CALL mpi_allgather(buffer, size1*loc_size, MPI_DOUBLE_PRECISION, dat, size1*loc_size, MPI_DOUBLE_PRECISION, mpi_comm_world, ierr)
  !
  deallocate(buffer)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx1(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: dat_size
  complex(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
  integer loc_size
  complex(dp), allocatable :: buffer(:)
  !
#if defined __MPI
  loc_size=dat_size/nnode
  allocate(buffer(1:loc_size))
  !
  buffer(:)=dat(inode*loc_size+1:(inode+1)*loc_size)
  !
  CALL mpi_allgather(buffer, loc_size, MPI_DOUBLE_COMPLEX, dat, loc_size, MPI_DOUBLE_COMPLEX, mpi_comm_world, ierr)
  !
  deallocate(buffer)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx3(dat, size1, size2, size3)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: size1, size2, size3
  complex(dp) :: dat(1:size1, 1:size2, 1:size3)
  !
  integer ierr, dat_size
  !
  integer loc_size
  complex(dp), allocatable :: buffer(:, :, :)
  !
#if defined __MPI
  dat_size=size1*size2*size3
  loc_size=size3/nnode
  allocate(buffer(1:size1, 1:size2, 1:loc_size))
  !
  buffer(:, :, :)=dat(:, :, inode*loc_size+1:(inode+1)*loc_size)
  !
  CALL mpi_allgather(buffer, size1*size2*loc_size, MPI_DOUBLE_COMPLEX, dat, size1*size2*loc_size, MPI_DOUBLE_COMPLEX, mpi_comm_world, ierr)
  !
  deallocate(buffer)
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
