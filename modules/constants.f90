MODULE constants
  !
  implicit none
  !
  integer, parameter :: dp=selected_real_kind(14, 200)
  real(dp), parameter :: twopi=3.141592653589793*2.d0
  real(dp), parameter :: sqrtpi=1.7724538509055158819194275566d0
  real(dp), parameter :: sqrt2=1.41421356237309504880168872421d0
  real(dp), parameter :: logpi_2=0.572364942924700087071713675675d0
  complex(dp), parameter :: cmplx_1=cmplx(1.d0, 0.d0, KIND=dp)
  complex(dp), parameter :: cmplx_i=cmplx(0.d0, 1.d0, KIND=dp)
  complex(dp), parameter :: cmplx_0=cmplx(0.d0, 0.d0, KIND=dp)
  integer, parameter :: stdin=5
  integer, parameter :: stdout=6
  integer, parameter :: fin=10
  integer, parameter :: fout=11
  integer, parameter :: maxint=2147483647
  real(dp), parameter :: eps4=1.0d-4
  real(dp), parameter :: eps6=1.0d-6
  real(dp), parameter :: eps9=1.0d-9
  !
END MODULE
