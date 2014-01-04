MODULE orbitals
  !
  use constants
  use wanndata, only: norb
  !
  IMPLICIT NONE
  !
  ! crystal structure related info
  !
  integer nsp, nat    ! # of species & total # of atoms
  integer, allocatable :: sp_nat(:)   ! # of atoms per specie
  integer, allocatable :: at_sp(:)   ! specie index of the atom
  real(dp), allocatable :: tau(:, :) ! atomic positions (fractional coordinates)
  !
  ! orbital related info
  !
  integer, allocatable :: gorb_at(:) ! global orbital index to atom index
  integer, allocatable :: gorb_loc(:) ! global orbital index to local (atomic) orbital index
  integer, allocatable :: glob_orb(:, :) ! global orbital index from atom index & local orbital index
  integer, allocatable :: sp_norb(:)  ! # of orbitals per specie
  !
CONTAINS
  !
SUBROUTINE read_orbitals
  !
  use constants
  use wanndata, only: norb
  !
  IMPLICIT NONE
  !
  integer ii, jj, kk, iatm, iorb
  !
  open(unit=fin, file="orbdef.in")
  !
  read(fin, *) nsp
  !
  allocate(sp_nat(1:nsp))
  allocate(sp_norb(1:nsp))
  !
  read(fin, *) sp_nat(:)
  read(fin, *) sp_norb(:)
  !
  if (norb.ne.SUM(sp_nat(:)*sp_norb(:))) then
    write(stdout, *) " !!! ERROR: inconsistent orbital definition and hamiltonian file."
    stop
  endif
  !
  nat=SUM(sp_nat(:))
  !
  write(stdout, *) " Reading orbital definitions of ", nat, " atoms"
  !
  allocate(at_sp(1:nat))
  allocate(tau(1:nat, 1:3))
  allocate(gorb_at(1:norb))
  allocate(gorb_loc(1:norb))
  allocate(glob_orb(1:nat, 1:maxval(sp_norb)))
  iorb=0
  iatm=0
  do ii=1, nsp
    do jj=1, sp_nat(ii)
      iatm=iatm+1
      at_sp(iatm)=ii
      read(fin, *) tau(iatm,1:3)
      do kk=1, sp_norb(ii)
        iorb=iorb+1
        gorb_at(iorb)=iatm
        gorb_loc(iorb)=kk
        glob_orb(iatm, kk)=iorb
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  !
END SUBROUTINE
  !
SUBROUTINE finalize_orbitals
  !
  IMPLICIT NONE
  !
  if(allocated(sp_nat)) deallocate(sp_nat)
  if(allocated(sp_norb)) deallocate(sp_norb)
  if(allocated(at_sp)) deallocate(at_sp)
  if(allocated(tau)) deallocate(tau)
  if(allocated(gorb_at)) deallocate(gorb_at)
  if(allocated(gorb_loc)) deallocate(gorb_loc)
  if(allocated(glob_orb)) deallocate(glob_orb)
  !
END SUBROUTINE
  !
END MODULE
