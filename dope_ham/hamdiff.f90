PROGRAM ham_diff
  !
  use constants
  use wanndata
  !
  IMPLICIT NONE
  !
  character*80 fn
  integer irpt, iorb, jorb
  !
  first_ham=.true.
  fn="ham0.dat"
  CALL read_ham(fn)
  fn="ham1.dat"
  CALL read_ham(fn)
  !
  open(unit=fout, file="hamdiff.dat")
  write(fout, *) "Hamiltonian difference"
  write(fout, *) norb
  write(fout, *) nrpt
  write(fout, '(15I5)') (weight(irpt), irpt=1, nrpt)
  do irpt=1,nrpt
    do iorb=1,norb
      do jorb=1,norb
        write(fout,'(5I5,2F12.6)') rvec(:, irpt), jorb, iorb, ham(1,iorb,jorb,irpt)-ham(0,iorb,jorb,irpt)
      enddo
    enddo
  enddo
  !
  close(unit=fout)

  !
  CALL finalize_wann
  !
END PROGRAM
