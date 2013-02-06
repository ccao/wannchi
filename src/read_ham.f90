SUBROUTINE read_ham(seed)
!
  USE constants
  USE wanndata
  !
  IMPLICIT NONE
  !
  CHARACTER(len=80) seed
  INTEGER irpt, iorb, jorb, t1, t2, t3, t4, t5
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  write(stdout, *) "Reading file "//trim(seed)//"_hr.dat"
  !
  open(unit=fin, file=trim(seed)//"_hr.dat")
  !
  read(fin, *)
  read(fin, *) norb
  read(fin, *) nrpt
  !
  write(stdout, *) "  Dimensions:"
  write(stdout, *) "    # of orbitals:", norb
  write(stdout, *) "    # of real-space grid:", nrpt
  !
  allocate(ham(1:norb, 1:norb, 1:nrpt))
  allocate(weight(1:nrpt))
  allocate(rvec(1:nrpt, 1:3))
  !
  allocate(wt(1:nrpt))
  read(fin, '(15I5)') (wt(irpt),irpt=1,nrpt)
  weight(:)=wt(:)
  deallocate(wt)
  !
  do irpt=1, nrpt
    do iorb=1, norb
      do jorb=1, norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if ((iorb.eq.1).and.(jorb.eq.1)) then
          rvec(irpt,1)=t1
          rvec(irpt,2)=t2
          rvec(irpt,3)=t3
        endif
        ham(iorb, jorb, irpt)=CMPLX(a,b)
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  write(stdout, *) "Done."
  !
END SUBROUTINE
