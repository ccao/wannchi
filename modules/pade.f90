MODULE pade_sum
  !
  use constants, only : dp, twopi
  !
  implicit none
    !
    integer npole
    real(dp), dimension(:), allocatable :: zp, eta
    !
  contains
  !
  SUBROUTINE init_pade(np_, ispade)
    !
    implicit none
     !
    integer np_
    logical ispade
    !
    real(dp), dimension(2*np_, 2*np_) :: A, B, V, BV, D, BA
    real(dp), dimension(18*np_) :: work
    real(dp), dimension(2*np_)  :: zp_, eta_
    !
    integer, dimension(2*np_) :: ipiv
    !
    integer ii
    !
    npole=np_
    allocate(zp(npole), eta(npole))
    !
    if (.not.ispade) then
      do ii=1, npole
        zp(ii)=(ii-0.5d0)*twopi
        eta(ii)=1.d0
      enddo
    else
      !
      B=0.d0
      D=0.d0
      !
      do ii=1, 2*np_
        D(ii, ii)=sqrt(2.d0*ii-1.d0)
      enddo
      !
      do ii=1, 2*np_-1
        B(ii, ii+1)=-0.5d0
        B(ii+1, ii)=-0.5d0
      enddo
      ! inv (B)
      call dgetrf(2*np_, 2*np_, B, 2*np_, ipiv, ii)
      call dgetri(2*np_, B, 2*np_, ipiv, work, 18*np_, ii)
      ! BV=B*D  BA=D*BV
      call dgemm('N', 'N', 2*np_, 2*np_, 2*np_, 1.d0, B, 2*np_, D, 2*np_, 0.d0, BV, 2*np_)
      call dgemm('N', 'N', 2*np_, 2*np_, 2*np_, 1.d0, D, 2*np_, BV, 2*np_, 0.d0, BA, 2*np_)
      ! zp_=eig(BA)
      call dsyev('V', 'U', 2*np_, BA, 2*np_, zp_, work, 18*np_, ii)
      !
      do ii=1, 2*np_
        D(ii, ii)=1.d0/sqrt(2.d0*ii-1.d0)
      enddo
      !  V=D*BA
      call dgemm('N', 'N', 2*np_, 2*np_, 2*np_, 1.d0, D, 2*np_, BA, 2*np_, 0.d0, V, 2*np_)
      !  eig(A)
      A(:, :) = V(:, :)
      call dgetrf(2*np_, 2*np_, A, 2*np_, ipiv, ii)
      call dgetri(2*np_, A, 2*np_, ipiv, work, 18*np_, ii)
      !  BV=A*B
      call dgemm('N', 'N', 2*np_, 2*np_, 2*np_, 1.d0, A, 2*np_, B, 2*np_, 0.d0, BV, 2*np_)
      ! V*BV
      do ii=1, 2*np_
        eta_(ii)=V(1, ii)*BV(ii, 1)*zp_(ii)/4.d0
      enddo
      !
      do ii=1, npole
        zp(ii)=(zp_(npole+ii)-zp_(npole+1-ii))/2.d0
        eta(ii)=(eta_(npole+ii)+eta_(npole+1-ii))/2.d0
      enddo
      !
    endif
    !
  END SUBROUTINE

  SUBROUTINE print_pade
    !
    use constants,  only : stdout, twopi
    use para,       only : inode
    !
    implicit none
    !
    integer ii
    !
    if (inode.eq.0) then
      write(stdout, '(A,1I5,A)')  "  #  Pade summation initialized with ", npole, " frequencies."
      write(stdout, '(A)') "  #   ii         Zp(ii)         Eta(ii)"
      do ii=1, npole
        write(stdout, '(A,1I5,2G16.9)') "  # ", ii, zp(ii)*2.d0/twopi, eta(ii)
      enddo
    endif
    !
  END SUBROUTINE

  SUBROUTINE finalize_pade
    !
    implicit none
    !
    if (allocated(zp)) deallocate(zp)
    if (allocated(eta)) deallocate(eta)
    !
  END SUBROUTINE
  !
END MODULE
  
