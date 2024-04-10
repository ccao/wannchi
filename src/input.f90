MODULE input
  !
  use constants
  !
  implicit none
  !
  character(len=80) seed
  ! SeedName
  !
  real(dp) mu, beta
  ! Fermi level and System temperature
  !
  integer nqpt
  ! Number of qpts
  !
  real(dp), dimension(:, :), allocatable :: qvec
  ! qvectors: qvec(3, iq)
  !
  integer nnu
  ! Number of frequencies to be calculated
  !
  real(dp) emin, emax
  ! (In case of spectrum calculation)
  ! emin, emax and nnu determines nu
  !
  real(dp) eps
  ! Input infinitesmal
  !
  complex(dp), dimension(:), allocatable :: nu
  ! In case of Matsubara calculation
  ! nu is determined by beta and nnu
  !
  logical spectra_calc
  ! Calculate on real axis? (spectra calculation)
  !
  logical trace_only
  ! Calculate only trace of Chi
  !
  logical ff_only
  ! Calculate only interaction part
  !
  logical use_lehman
  ! Use Lehman Representation To Calculate
  !
  logical fast_calc
  ! Use fast algorithm (more memory required)
  !
  integer npade
  ! Number of Poles in Pade Summation
  !
  namelist /CONTROL/ use_lehman, trace_only, ff_only, fast_calc, eps, nnu, emin, emax, npade
  namelist /SYSTEM/ spectra_calc, seed, beta, mu
  !
CONTAINS
  !
 SUBROUTINE read_input(codename)
  !
  !
  use constants, only : dp, eps6, fin, fout
  use para,      only : inode, para_sync_int, para_sync_real
  !
  implicit none
  !
  character(*) codename
  !
  character(len=80) line
  !
  integer, dimension(8)  :: tt_int
  real(dp), dimension(5) :: tt_real
  integer ii
  !
  ! Example INPUT FILE:
  ! &SYSTEM
  !   seed='wannier90',
  !   beta=2000.d0,
  !   mu=0.d0,
  !   spectra_calc = .false.
  ! /
  ! &CONTROL
  !   use_lehman=.false.
  !   trace_only=.false.
  !   ff_only  = .true.
  !   fast_calc=.true.
  !   npade=80
  !   nnu=1
  !   ! eps=0.001  ! Only Used for Spectra
  !   ! emin=0.0   !
  !   ! emax=0.0
  ! /
  !
  ! DEFAULT:
  !
  seed='wannier90'
  ! Default seed name
  !
  beta=1e7
  ! Very low T calculation (~0K)
  !
  mu=0.d0
  ! Default Fermi level 
  !
  spectra_calc = .false.
  ! Calculate on imaginary frequency
  !
  use_lehman = .false.
  ! Use G*G algorithm
  !
  trace_only = .false.
  ! Only calculate trace
  !
  ff_only    = .true.
  ! Calculate only the interaction part
  !
  fast_calc = .true.
  ! Use more memory
  !
  npade=80
  ! 80 Pade Poles
  !
  nnu=1
  ! Calculate only single frequency
  !
  eps=eps6
  ! Default small imaginary part
  !
  emin=0.0
  ! Default at Ef
  !
  emax=0.0
  ! 
  ! logical spectra_calc, use_lehman, trace_only, ff_only, fast_calc
  ! real    eps, emin, emax, beta, mu
  ! integer nnu
  !
  if (inode.eq.0) then
    ! Read Structure Input
    !
    open(unit=fin, file=trim(codename)//".inp")
    !
    read(nml=SYSTEM, unit=fin)
    read(nml=CONTROL, unit=fin)
    !
    close(unit=fin)
    !
    tt_int(:)=0
    !
    if (spectra_calc) tt_int(1) = 1
    if (use_lehman)   tt_int(2) = 1
    if (trace_only)   tt_int(3) = 1
    if (ff_only)      tt_int(4) = 1
    if (fast_calc)    tt_int(5) = 1
    tt_int(6) = nnu
    tt_int(7) = npade
    !
    tt_real(:)=(/eps, emin, emax, beta, mu/);
    !
  endif
  !
  call para_sync_int(tt_int, 7)
  !
  spectra_calc = (tt_int(1).eq.1)
  use_lehman   = (tt_int(2).eq.1)
  trace_only   = (tt_int(3).eq.1)
  ff_only      = (tt_int(4).eq.1)
  fast_calc    = (tt_int(5).eq.1)
  nnu          =  tt_int(6)
  npade        =  tt_int(7)
  !
  call para_sync_real(tt_real, 5)
  eps  = tt_real(1)
  emin = tt_real(2)
  emax = tt_real(3)
  beta = tt_real(4)
  mu   = tt_real(5)
  !
  allocate(nu(nnu))
  !
  if (spectra_calc) then
    !
    do ii=1, nnu
      nu(ii)=emin+(emax-emin)*(ii-1)/(nnu-1)
    enddo
    !
  else
    !
    do ii=1, nnu
      nu(ii)=twopi*(ii-1)/beta
    enddo
    !
  endif
  !
  if (beta<0) beta=1.d7
  if (eps>eps4.or.eps<eps9) eps=eps6
  !
 END SUBROUTINE
  !
 SUBROUTINE read_qpoints()
  !
  use constants, only : dp, fin
  use para,      only : inode, para_sync_int, para_sync_real
  !
  implicit none
  !
  integer mode
  integer nq1, nq2, nq3
  integer iq1, iq2, iq3
  real(dp), dimension(3) :: tq1, tq2, tq0
  integer, dimension(4) :: tt
  !
  if (inode.eq.0) then
    !
    open(unit=fin, file="QPOINTS")
    !
    ! QPOINTS Example 1
    !   0           ! Single Point Calculation
    !  0.5 0.5 0.5  ! Qvec
    !
    ! QPOINTS Example 2
    !   1           ! Line mode
    !   1  48       !  nseg   ninterpolate
    ! 0.0  0.0  0.0    0.5  0.0  0.0  ! seg 1: Q1  Q2
    !
    ! QPOINTS Example 3
    !   2           ! Plane mode
    !  48  48       ! nint1  nint2
    ! 0.0  0.0  0.0 ! Vertex
    ! 1.0  0.0  0.0 ! Direction 1
    ! 0.0  1.0  0.0 ! Direction 2
    !
    ! QPOINTS Example 4
    !   3           ! Full BZ
    !  48  48  48   ! nint1  nint2  nint3
    !
    read(fin, *) mode
    !
    if (mode.eq.0) then
      nqpt=1
      nq1=1
      nq2=1
      nq3=1
    elseif (mode.eq.1) then
      nq3=1
      read(fin, *) nq1, nq2
      nqpt=nq1*(nq2+1)
    elseif (mode.eq.2) then
      nq3=1
      read(fin, *) nq1, nq2
      nqpt=(nq1+1)*(nq2+1)
    elseif (mode.eq.3) then
      read(fin, *) nq1, nq2, nq3
      nqpt=nq1*nq2*nq3
    endif
    !
    tt(1)=nq1
    tt(2)=nq2
    tt(3)=nq3
    tt(4)=nqpt 
    !
  endif
  !
  call para_sync_int(tt, 4)
  !
  nq1=tt(1)
  nq2=tt(2)
  nq3=tt(3)
  nqpt=tt(4)
  !
  allocate(qvec(3, nqpt))
  !
  if (inode.eq.0) then
    !
    if (mode.eq.0) then
      read(fin, *) qvec(:, 1)
    elseif (mode.eq.1) then
      do iq1=0, nq1-1
        read(fin, *) tq1(:), tq2(:)
        do iq2=0, nq2
          qvec(:, iq1*(nq2+1)+iq2+1)=tq1 + (iq2*1.d0)/nq2*tq2
        enddo
      enddo
    elseif (mode.eq.2) then
      read(fin, *) tq0
      read(fin, *) tq1
      read(fin, *) tq2
      !
      do iq1=0, nq1
        do iq2=0, nq2
          qvec(:, iq1*(nq2+1)+iq2+1)=tq0 + (iq1*1.d0)/nq1*tq1 + (iq2*1.d0)/nq2*tq2
        enddo
      enddo
    elseif (mode.eq.3) then
      do iq1=0, nq1-1
        do iq2=0, nq2-1
          do iq3=0, nq3-1
            qvec(1, iq1*nq2*nq3+iq2*nq3+iq3+1) = iq1*1.d0/nq1
            qvec(2, iq1*nq2*nq3+iq2*nq3+iq3+1) = iq2*1.d0/nq2
            qvec(3, iq1*nq2*nq3+iq2*nq3+iq3+1) = iq3*1.d0/nq3
          enddo
        enddo
      enddo
    endif
    !
    close(fin)
    !
  endif
  !
  call para_sync_real(qvec, nqpt*3)
  !
 END SUBROUTINE

  !
 SUBROUTINE finalize_input
  !
  implicit none
  !
  if (allocated(nu))   deallocate(nu)
  if (allocated(qvec)) deallocate(qvec)
  !
 END SUBROUTINE
  !
END MODULE
