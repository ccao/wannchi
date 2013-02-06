MODULE constants
  !
  implicit none
  !
  integer, parameter :: dp=selected_real_kind(14, 200)
  real(dp), parameter :: twopi=3.141592653589793*2.d0
  complex(dp), parameter :: cmplx_i=cmplx(0.d0, 1.d0)
  complex(dp), parameter :: cmplx_0=cmplx(0.d0, 0.d0)
  integer, parameter :: stdout=6
  integer, parameter :: fin=10
  integer, parameter :: fout=11
  real(dp), parameter :: eps=1.0d-6
  !
END MODULE


