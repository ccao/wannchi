PROGRAM mixham
!
  USE constants
  !
  IMPLICIT NONE
  !
  CHARACTER(len=80) fn
  !
  INTEGER irpt, iorb, jorb, t1, t2, t3, t4, t5
  REAL(dp) a, b
  REAL(dp) wt
  INTEGER, ALLOCATABLE :: weight(:)
  INTEGER, ALLOCATABLE :: rvec(:, :)
  INTEGER nrpt, norb
  COMPLEX(dp), ALLOCATABLE :: ham(:,:,:)
  !
  write(stdout, *) " Weight of first hamiltonian:"
  read(stdin, *) wt
  write(stdout, *) "  First hamiltonian file:"
  read(stdin, *) fn
  !
  open(unit=fin, file=fn)
  !
  read(fin, *)
  read(fin, *) norb
  read(fin, *) nrpt
  !
  allocate(ham(1:norb, 1:norb, 1:nrpt))
  allocate(weight(1:nrpt))
  allocate(rvec(1:3, 1:nrpt))
  !
  read(fin, '(15I5)') (weight(irpt), irpt=1, nrpt)
  do irpt=1, nrpt
    do iorb=1, norb
      do jorb=1, norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if ((iorb.eq.1).and.(jorb.eq.1)) then
          rvec(1, irpt)=t1
          rvec(2, irpt)=t2
          rvec(3, irpt)=t3
        endif
        ham(iorb, jorb, irpt)=CMPLX(a, b)*wt
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  !
  write(stdout, *) "  Second hamiltonian file:"
  read(stdin, *) fn
  !
  open(unit=fin, file=fn)
  !
  read(fin, *)
  read(fin, *) norb
  read(fin, *) nrpt
  !
  read(fin, '(15I5)') (weight(irpt), irpt=1, nrpt)
  do irpt=1, nrpt
    do iorb=1, norb
      do jorb=1, norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if ((iorb.eq.1).and.(jorb.eq.1)) then
          rvec(1, irpt)=t1
          rvec(2, irpt)=t2
          rvec(3, irpt)=t3
        endif
        ham(iorb, jorb, irpt)=ham(iorb, jorb, irpt)+(1.d0-wt)*CMPLX(a, b)
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  !
  open(unit=fout, file="newham.dat")
  !
  write(fout, *) "Mixed Hamiltonian:"
  write(fout, *) norb
  write(fout, *) nrpt
  write(fout, '(15I5)') (weight(irpt), irpt=1, nrpt)
  do irpt=1,nrpt
    do iorb=1,norb
      do jorb=1,norb
        write(fout,'(5I5,2F12.6)') rvec(:,irpt), jorb, iorb, ham(jorb,iorb,irpt)
      enddo
    enddo
  enddo
  !
  close(unit=fout)
  !
  deallocate(ham)
  deallocate(weight)
  deallocate(rvec)
END PROGRAM
