!
!   wanndata.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!

MODULE wanndata
  !
  ! This module encapsules the operations on Wannier Hamiltonian
  !   BRIEFS:
  !     Initialize 
  !         by calling   read_ham
  !       Dimension of _hr.dat can be obtained 
  !         by calling   read_ham_dim
  !     Finalize 
  !         by calling   finalize_wann
  !     (Irreversible) Fermi level shift 
  !         by calling   wannham_shift_ef
  !     Obtain K-space hk
  !         by calling   calc_hk
  !
  !
  use constants, only : dp
  !
  IMPLICIT NONE
  !
  TYPE wannham
    !
    INTEGER norb
    !
    REAL(DP), ALLOCATABLE :: tau(:,:)
    !
    INTEGER nrpt
    INTEGER r000 ! This is the on-site
    !
    COMPLEX(DP), ALLOCATABLE :: hr(:,:,:)
    !
    REAL(DP), ALLOCATABLE :: weight(:)
    !
    REAL(DP), ALLOCATABLE :: rvec(:,:)
    !
  END TYPE ! wannham
  !
CONTAINS

SUBROUTINE wannham_shift_ef(ham, mu)
  !
  USE constants,  only : dp
  !
  implicit none
  !
  real(dp), intent(in)        :: mu
  TYPE(wannham)               :: ham
  !
  integer ii
  !
  do ii=1, ham%norb
    ham%hr(ii, ii, ham%r000) = ham%hr(ii, ii, ham%r000)-mu
  enddo
  !
END SUBROUTINE

SUBROUTINE write_ham(ham, seed)
  !
  USE constants, only : stdout, fout
  USE para,      only : inode
  !
  IMPLICIT NONE
  !
  CHARACTER(*), intent(in)   :: seed
  TYPE(wannham), intent(in)  :: ham
  !
  INTEGER irpt, iorb, jorb
  !
  if (inode.eq.0) then
    write(stdout, *) " # Write patched Wannier Hamiltonian "//trim(seed)//"_hr.dat"
    !
    open(unit=fout, file=trim(seed)//"_hr.dat")
    !
    write(fout, *) "# Patched Wannier Hamiltonian"
    write(fout, '(1I10)') ham%norb
    write(fout, '(1I10)') ham%nrpt
    !
    write(fout, '(15I5)') nint(ham%weight(:))
    do irpt=1, ham%nrpt
      do iorb=1, ham%norb
        do jorb=1, ham%norb
          write(fout, '(5I5,2F22.16)') nint(ham%rvec(:, irpt)), jorb, iorb, ham%hr(jorb, iorb, irpt)
        enddo
      enddo
    enddo
    !
    close(unit=fout)
    write(stdout, *) "Done."
  endif
  !
END SUBROUTINE

SUBROUTINE read_ham_dim(ham, seed)
  !
  use para,      only : para_sync_int, inode
  use constants, only : fin, stdout
  !
  implicit none
  !
  CHARACTER(*), intent(in)   :: seed
  TYPE(wannham), intent(out) :: ham
  !
  integer,dimension(5) :: tt
  if (inode.eq.0) then
    write(stdout, *) " # Hamiltonian dimension read from file "//trim(seed)//"_hr.dat"
    !
    open(unit=fin, file=trim(seed)//"_hr.dat")
    !
    read(fin, *)
    read(fin, *) tt(1) ! norb
    read(fin, *) tt(2) ! nrpt
    !
    close(unit=fin)
    !
  endif
  !
  CALL para_sync_int(tt, 2)
  ham%norb=tt(1)
  ham%nrpt=tt(2)
  allocate(ham%tau(3, ham%norb))
  !
  ham%tau(:,:)=0.d0
  !
  if (inode.eq.0) then
    write(stdout, *) " #  Dimensions:"
    write(stdout, *) "    # of orbitals:", ham%norb
    write(stdout, *) "    # of real-space grid:", ham%nrpt
  endif
  !
END SUBROUTINE

SUBROUTINE read_ham(ham, seed)
  !
  USE para,      only: para_sync_int, para_sync_real, para_sync_cmplx, inode, para_sync0
  USE constants, only: dp, fin, stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(*), intent(in)   :: seed
  TYPE(wannham)              :: ham
  !
  INTEGER irpt, iorb, jorb
  integer, dimension(5) :: tt
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  if (inode.eq.0) then
    write(stdout, *) " # Reading file "//trim(seed)//"_hr.dat"
    !
    open(unit=fin, file=trim(seed)//"_hr.dat")
    !
    read(fin, *)
    read(fin, *) tt(1) ! norb
    read(fin, *) tt(2) ! nrpt
    !
    write(stdout, *) " #  Dimensions:"
  endif
  !
  CALL para_sync_int(tt, 2)
  ham%norb=tt(1)
  ham%nrpt=tt(2)
  !
  allocate(ham%hr(ham%norb, ham%norb, ham%nrpt))
  allocate(ham%weight(ham%nrpt))
  allocate(ham%rvec(3, ham%nrpt))
  allocate(ham%tau(3, ham%norb))
  !
  ! site position initialized to 0
  !   They are to be read from POSCAR file
  !
  ham%tau(:,:)=0.d0
  !
  if (inode.eq.0) then
    write(stdout, *) "    # of orbitals:", ham%norb
    write(stdout, *) "    # of real-space grid:", ham%nrpt
    allocate(wt(1:ham%nrpt))
    read(fin, '(15I5)') (wt(irpt),irpt=1,ham%nrpt)
    ham%weight(:)=wt(:)
    deallocate(wt)
    !
    do irpt=1, ham%nrpt
      do iorb=1, ham%norb
        do jorb=1, ham%norb
          read(fin, *) tt, a, b
          if ((jorb.eq.1).and.(iorb.eq.1)) then
            ham%rvec(:, irpt)=tt(1:3)
            if (tt(1)**2+tt(2)**2+tt(3)**2.eq.0) then
              ham%r000=irpt
            endif
          endif
          ham%hr(jorb, iorb, irpt)=CMPLX(a, b, KIND=dp)
        enddo
      enddo
    enddo
    !
    close(unit=fin)
    write(stdout, *) " # Done."
  endif
  !
  CALL para_sync_cmplx(ham%hr, ham%norb * ham%norb * ham%nrpt)
  CALL para_sync_real(ham%weight, ham%nrpt)
  CALL para_sync_real(ham%rvec, 3*ham%nrpt)
  CALL para_sync0(ham%r000)
  !
  !if(inode.eq.0) then
  !  !
  !  open(unit=fin, file="tau.xyz")
  !  read(fin, *) tt(1), tt(3)
  !  jorb=0
  !  do iorb=1, tt(1)
  !    read(fin, *) tt(2), xx
  !    do irpt=1, tt(2)
  !      tau(:, jorb+1)=xx(:)
  !      jorb=jorb+1
  !    enddo
  !  enddo
  !  close(unit=fin)
  !  !
  !  if (tt(3).eq.1) then
  !    ! spinor is true
  !    do iorb=1, jorb
  !      tau(:, jorb+iorb)=tau(:, iorb)
  !    enddo
  !    jorb=jorb*2
  !  endif
  !  !
  !  if (jorb.ne.norb) then
  !    write(stdout, *) "!!! FATAL ERROR: # of orbitals doesn't match!", jorb, norb
  !    stop
  !  endif
  !  !
  !endif
  !!
  !CALL para_sync_real(tau, 3*norb)
  !!
  !do irpt=1, nrpt
  !  if (sum(rvec(:, irpt)**2)<eps4) r000=irpt
  !enddo
  !
END SUBROUTINE

SUBROUTINE finalize_wann(ham, all)
  !
  IMPLICIT NONE
  !
  TYPE(wannham), intent(inout) :: ham
  logical, intent(in)          :: all
  !
  if (allocated(ham%hr)) deallocate(ham%hr)
  if (all.and.allocated(ham%weight)) deallocate(ham%weight)
  if (all.and.allocated(ham%rvec)) deallocate(ham%rvec)
  if (all.and.allocated(ham%tau)) deallocate(ham%tau)
  !
END SUBROUTINE

SUBROUTINE calc_hk(hk, ham, kvec)
  !
  use constants,     only: dp, twopi, cmplx_0
  !
  implicit none
  !
  TYPE(wannham), intent(in)          :: ham
  real(dp), dimension(3), intent(in) :: kvec
  complex(dp), dimension(ham%norb, ham%norb), intent(out) :: hk
  !logical, optional :: tauphase_correct
  !
  integer ir
  real(dp)    :: rdotk
  complex(dp) :: fact
  !
  integer io, jo
  real(dp)    :: ktau
  complex(dp) :: orbfac
  complex(dp), dimension(ham%norb)  :: phase
  !
  !logical additional_phase
  !
  hk(:,:)=cmplx_0
  !
  !if (.not. PRESENT(tauphase_correct)) then
  !  additional_phase=.true.
  !else
  !  additional_phase=tauphase_correct
  !endif
  !
  !if (additional_phase) then
  do io=1, ham%norb
    ktau=sum(kvec(:)*ham%tau(:, io))*twopi
    phase(io)=cmplx(cos(ktau), sin(ktau), KIND=dp)
  enddo
  !endif
  !
  do ir=1, ham%nrpt
    rdotk=sum(kvec(:)*ham%rvec(:, ir))*twopi
    fact=cmplx(cos(rdotk), sin(rdotk), KIND=dp)/ham%weight(ir)
    !
    do io=1, ham%norb
      do jo=1, ham%norb
        !
        !if (additional_phase) then
        hk(io, jo)=hk(io, jo)+fact*conjg(phase(io))*phase(jo)*ham%hr(io, jo, ir)
        !else
        !  hk(io, jo)=hk(io, jo)+fact*ham%hr(io, jo, ir)
        !endif
        !
      enddo
    enddo
  enddo
  !
END SUBROUTINE

END MODULE

