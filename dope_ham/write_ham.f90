SUBROUTINE write_ham
  !
  use constants, only: fout
  use wanndata, only: norb, nrpt, rvec, weight
  use dopedata, only: doped_ham
  !
  IMPLICIT NONE
  !
  integer irpt, iorb, jorb
  !
  open(unit=fout, file="newham.dat")
  write(fout, *) "New Hamiltonian"
  write(fout, *) norb
  write(fout, *) nrpt
  write(fout, '(15I5)') (weight(irpt), irpt=1, nrpt)
  do irpt=1,nrpt
    do iorb=1,norb
      do jorb=1,norb
        write(fout,'(5I5,2F12.6)') rvec(:, irpt), jorb, iorb, doped_ham(iorb,jorb,irpt)
      enddo
    enddo
  enddo
  !
  close(unit=fout)
  !
END SUBROUTINE
