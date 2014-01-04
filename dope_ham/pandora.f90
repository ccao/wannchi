PROGRAM pandora
  !
  use constants
  use dopedata
  use orbitals
  use wanndata
  !
  IMPLICIT NONE
  !
  character*80 fn
  integer fid
  !
  first_ham=.true.
  fid=16
  open(unit=fid, file="pandora.in")
  read(fid, *) fn
  CALL read_ham(fn)
  read(fid, *) fn
  CALL read_ham(fn)
  !
  CALL read_orbitals
  !
  read(fid, *) ndopant
  read(fid, *) ref_dopant
  !
  CALL initialize_dopedata
  !
  read(fid, *) dopant_list(:)
  !
  write(stdout, *) " ... Will perform hamiltonian translation on ..."
  write(stdout, *) "      ", ndopant, " atomic sites"
  write(stdout, *) dopant_list
  !
  write(stdout, *) " !!!: Reference dopant coordinate:"
  write(stdout, *) tau(ref_dopant, :)
  !
  close(unit=fid)
  !
  write(stdout, *) " ... Translating and mixing hamiltonian ..."
  !
  CALL mix_ham
  !
  CALL write_ham
  !
  CALL finalize_dopedata
  CALL finalize_orbitals
  CALL finalize_wann
  !
END PROGRAM
