REAL FUNCTION calc_occ(en)
  !
  use constants, only : dp
  use banddata, only : ef
  use input, only : temp
  !
  implicit none
  !
  if ( temp < eps6 ) then ! Zero temperature
    if ( en > ef ) then
      calc_occ=0.d0
    else if ( en < ef ) then
      calc_occ=1.d0
    else
      calc_occ=0.5d0
    endif
  else  ! Fermi-Dirac distribution for non-zero temperature
    calc_occ=1.d0/(1.d0+exp((en-ef)/temp))
  endif
  !
END FUNCTION
