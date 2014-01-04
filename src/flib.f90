<<<<<<< HEAD
FUNCTION fermi_dirac(en, temp)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) en, temp, fermi_dirac
  !
  fermi_dirac=1.d0/(1.d0+exp(en/temp))
  !
END FUNCTION

FUNCTION converged(diff, n)
  !
  use constants, only: dp, eps
  !
  integer, intent(in) :: n
  complex(dp), intent(in) :: diff(n,n,n,n)
  logical converged
  !
  integer i, j, k, l
  !
  converged = .true.
  do i=1, n
    do j=1, n
      do k=1, n
        do l=1, n
          if ( abs(diff(i, j, k, l)) < 1.0d-5 ) then
            converged = .false.
          endif
        enddo
      enddo
    enddo
  enddo
  !
END FUNCTION

SUBROUTINE matrix_multiply(c, a, b, n)
  !
  use constants, only: dp, cmplx_0
  !
  implicit none
  !
  integer, intent(in) :: n
  complex(dp), intent(out) :: c(n, n, n, n)
  complex(dp), intent(in) :: a(n, n, n, n), b(n, n, n, n)
  !
  integer ii, jj, kk, ll, il, jl
  !
  c = cmplx_0
  !
  do ii=1, n
    do jj=1, n
      do kk=1, n
        do ll=1, n
          do il=1, n
            do jl=1, n
              c(ii, jj, ll, kk) = c(ii, jj, kk, ll) + a(ii, jj, il, jl)*b(il, jl, kk, ll)
            enddo
          enddo
        enddo ! ll
      enddo
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE map_to_kpt(kv, ik, nkx, nky, nkz)
=======
INTEGER FUNCTION kpt_index(kv)
>>>>>>> New modulized version
  !
  use banddata, only : kvec, nkpt
  use constants, only : dp, eps6
  !
  implicit none
  !
  integer ik
  logical found
  real(dp) dk(1:3)
  !
  do ik=1, nkpt
    dk(:)=(kv(:)-kvec(:))-nint(kv(:)-kvec(:))
    if ( abs(dk(1))<eps6.and.abs(dk(2))<eps6.and.abs(dk(3))<eps6 ) then
      kpt_index=ik
      return
    endif
  enddo
  !
<<<<<<< HEAD
  ikx=mod(ik-1, nkx)
  iky=mod( (ik-1-ikx)/nkx, nky )
  ikz=(ik-1-iky*nkx-ikx)/(nkx*nky)
  kv(1)=ikx*(1.d0/nkx)
  kv(2)=iky*(1.d0/nky)
  kv(3)=ikz*(1.d0/nkz)
!  if (kv(1)>0.5) kv(1)=kv(1)-1.d0
!  if (kv(2)>0.5) kv(2)=kv(2)-1.d0
!  if (kv(3)>0.5) kv(3)=kv(3)-1.d0
=======
  kpt_index=-1
  return
>>>>>>> New modulized version
  !
END FUNCTION
