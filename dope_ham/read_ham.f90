SUBROUTINE read_ham(fn)
!
  USE constants
  USE wanndata
  !
  IMPLICIT NONE
  !
  CHARACTER*80 fn
  INTEGER irpt, iorb, jorb, t1, t2, t3, t4, t5
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  write(stdout, *) " # Reading file "//fn
  !
  open(unit=fin, file=fn)
  !
  read(fin, *)
  read(fin, *) t1
  read(fin, *) t2
  !
  if(first_ham) then
    norb=t1
    nrpt=t2
  else
    if((norb.ne.t1).or.(nrpt.ne.t2)) then
      write(stdout, *) " !!! ERROR: inconsistent hamiltonian dimensions"
      return
    endif
  endif
  !
  write(stdout, *) " #  Dimensions:"
  write(stdout, *) "    # of orbitals:", norb
  write(stdout, *) "    # of real-space grid:", nrpt
  !
  if(first_ham) then
    allocate(ham(0:1, 1:norb, 1:norb, 1:nrpt))
    allocate(weight(1:nrpt))
    allocate(rvec(1:3, 1:nrpt))
  endif
  !
  allocate(wt(1:nrpt))
  read(fin, '(15I5)') (wt(irpt),irpt=1,nrpt)
  if(first_ham) weight(:)=wt(:)
  deallocate(wt)
  !
  do irpt=1, nrpt
    do iorb=1, norb
      do jorb=1, norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if (first_ham.and.(iorb.eq.1).and.(jorb.eq.1)) then
          rvec(1, irpt)=t1
          rvec(2, irpt)=t2
          rvec(3, irpt)=t3
        endif
        if (first_ham) then
          ham(0, iorb, jorb, irpt)=CMPLX(a,b)
        else
          ham(1, iorb, jorb, irpt)=CMPLX(a,b)
        endif
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  write(stdout, *) " # Done."
  !
  if(first_ham) then
    first_ham=.false.
    !  else
    !    ham(1, :, :, :) = ham(1, :, :, :)-ham(0, :, :, :)
  endif
  !
END SUBROUTINE
