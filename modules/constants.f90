!
!   constants.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2023 CC. All rights reserved.
!
MODULE constants
  !
  implicit none
  !
  ! Data Kind Selector
  integer, parameter :: dp=selected_real_kind(14, 200)
  !
  ! Some useful constants
  real(dp), parameter :: twopi=6.283185307179586_dp
  real(dp), parameter :: sqrtpi=1.772453850905516_dp
  real(dp), parameter :: sqrt2=1.414213562373095_dp
  real(dp), parameter :: logpi_2=0.572364942924700_dp
  !
  ! Complex constants
  complex(dp), parameter :: cmplx_1=cmplx(1.d0, 0.d0, KIND=dp)
  complex(dp), parameter :: cmplx_i=cmplx(0.d0, 1.d0, KIND=dp)
  complex(dp), parameter :: cmplx_0=cmplx(0.d0, 0.d0, KIND=dp)
  !
  ! I/O units
  integer, parameter :: stdin=5
  integer, parameter :: stdout=6
  integer, parameter :: fin=10
  ! Primary input
  integer, parameter :: fout=11
  ! Primary output
  integer, parameter :: fout2=12
  ! Secondary output
  integer, parameter :: fin3=13
  integer, parameter :: fin4=14
  integer, parameter :: fin5=15
  integer, parameter :: fin6=16
  !
  integer, parameter :: fout3=17
  integer, parameter :: fout4=18
  integer, parameter :: fout5=19
  integer, parameter :: fout6=20
  integer, parameter :: fdebug=30
  ! Internal units...
  !
  ! Limits
  integer, parameter :: maxint=2147483647
  real(dp), parameter :: eps4=1.0d-4
  real(dp), parameter :: eps6=1.0d-6
  real(dp), parameter :: eps9=1.0d-9
  real(dp), parameter :: eps12=1.0d-12
  real(dp), parameter :: eps18=1.0d-18
  real(dp), parameter :: eps36=1.0d-36
  !
END MODULE
