!
! symmetry.f90
!
!   This file defines symmetry operations
!
  !   BRIEFS:
  !     Initialize 
  !         by calling   init_symm with angle/axis or matrix
  !     Generate <Ylm | C>
  !         by calling   generate_Ylm2C
  !     Generate <Ylm | L | Ylm'>
  !         by calling   generate_Lmatrix
  !     Generate <sigma | S | sigma'>
  !         by calling   generate_Smatrix
  !     Vector rotation
  !         by calling   rotate_vector
  !     Ylm basis rotation
  !         by calling   rotate_Ylm
  !     Ylms basis rotation
  !         by calling   rotate_Ylms
  !     Cubic basis rotation
  !         by calling   rotate_cubic
  !     Finalize not required 
  !
MODULE symmetry_module
  !
  use constants
  !
  implicit none
  !
  TYPE symmetry
    !
    real(dp), dimension(3, 3)   :: rot=reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/3,3/))
    ! Rotation matrix in Direct unit (lattice vector unit)
    real(dp), dimension(3)      :: tau=(/0.d0, 0.d0, 0.d0/)
    ! Translation matrix in Direct unit
    ! Matrix representation
    !
    real(dp), dimension(3)      :: axis=(/0.d0, 0.d0, 1.d0/)
    ! Axis in cartesian unit
    real(dp)                    :: theta=0.d0
    logical                     :: inv=.false.
    !
  END TYPE
  !
  interface init_symm
    module procedure init_symm_axis_angle_inv, init_symm_matrix
  end interface
  !
  CONTAINS
  !
  SUBROUTINE inverse_symm(symm)
    !
    use linalgwrap,  only: invmat
    !
    TYPE(symmetry), intent(inout)  :: symm
    !
    integer ii
    real(dp), dimension(3) :: tt
    !
    tt(:)=symm%tau(:)
    !
    CALL invmat(symm%rot, 3)
    !
    do ii=1, 3
      symm%tau(ii)=-sum(symm%rot(ii,:)*tt(:))
    enddo
    symm%theta=-symm%theta
    !
  END SUBROUTINE

  SUBROUTINE generate_Smatrix(Sx, Sy, Sz)
    !
    ! Pauli Matrices
    !
    complex(dp), dimension(2, 2), intent(out)   :: Sx, Sy, Sz
    ! Sx, Sy, Sz
    !
    Sx=cmplx_0
    Sy=cmplx_0
    Sz=cmplx_0
    !
    Sx(1,2)=cmplx_1
    Sx(2,1)=cmplx_1
    !
    Sy(1,2)=-cmplx_i
    Sy(2,1)=cmplx_i
    !
    Sz(1,1)=cmplx_1
    Sz(2,2)=-cmplx_1
    !
  END SUBROUTINE

  SUBROUTINE generate_Lmatrix(Lx, Ly, Lz, l)
    !
    ! <Ylm'|L|Ylm>
    !
    integer, intent(in)          :: l
    ! Angular momentum l
    complex(dp), dimension(2*l+1, 2*l+1), intent(out)   :: Lx, Ly, Lz
    ! Lx, Ly, Lz
    !
    complex(dp), dimension(2*l+1, 2*l+1) :: Lp, Lm
    !   Temporary matrices (L+, L-) :: Lp, Lm
    integer m
    !
    Lp=cmplx_0
    Lm=cmplx_0
    Lz=cmplx_0
    !
    do m=-l, l
      Lz(m+l+1, m+l+1)=m
      if (m<l) then
        Lp(m+l+2, m+l+1)=sqrt(real((l-m)*(l+m+1),dp))
        Lm(m+l+1, m+l+2)=sqrt(real((l+m+1)*(l-m),dp))
      endif
    enddo
    !
    Lx=(Lp+Lm)/2.d0
    Ly=(Lp-Lm)/(2.d0*cmplx_i)
    !
  END SUBROUTINE

  SUBROUTINE generate_Ylm2C(Umat, l)
    !
    ! This function generates <Ylm | C>
    !
    integer, intent(in)         :: l
    complex(dp), dimension(2*l+1, 2*l+1), intent(out)  :: Umat
    !
    integer  m
    !
    Umat=cmplx_0
    Umat(l+1, 1)=cmplx_1
    do m=1, l
      if (mod(m, 2).eq.1) then
        Umat(l+1-m, 2*m)=cmplx_1/sqrt2
        Umat(l+1+m, 2*m)=-cmplx_1/sqrt2
        Umat(l+1-m, 2*m+1)=cmplx_i/sqrt2
        Umat(l+1+m, 2*m+1)=cmplx_i/sqrt2
      else
        Umat(l+1-m, 2*m)=cmplx_1/sqrt2
        Umat(l+1+m, 2*m)=cmplx_1/sqrt2
        Umat(l+1-m, 2*m+1)=cmplx_i/sqrt2
        Umat(l+1+m, 2*m+1)=-cmplx_i/sqrt2
      endif
    enddo
    !
  END SUBROUTINE

  SUBROUTINE rotate_vector(v1, v0, symm, realspace)
    !
    real(dp), dimension(3), intent(out) :: v1
    real(dp), dimension(3), intent(in)  :: v0
    TYPE(symmetry), intent(in)          :: symm
    logical, optional                   :: realspace
    !
    integer ii
    !
    ! Default realspace to T
    if (.not.present(realspace)) then
      realspace=.true.
    endif
    !
    do ii=1, 3
      v1(ii)=sum(symm%rot(:, ii)*v0(:))
    enddo
    !
    if (realspace) then
      v1(:)=v1(:)+symm%tau(:)
    endif
    !
  END SUBROUTINE

  SUBROUTINE rotate_Ylms(rot, l, symm)
    !
    integer, intent(in) :: l
    TYPE(symmetry), intent(in) :: symm
    complex(dp), dimension(4*l+2, 4*l+2), intent(out) :: rot
    !
    complex(dp), dimension(2*l+1, 2*l+1) :: rotYlm
    complex(dp), dimension(2, 2) :: rotSpin
    !
    CALL rotate_spinor(rotSpin, symm)
    CALL rotate_Ylm(rotYlm, l, symm)
    !
    rot=cmplx_0
    rot(1:2*l+1, 1:2*l+1)=rotSpin(1,1)*rotYlm(:,:)
    rot(1:2*l+1, 2*l+2:4*l+2)=rotSpin(1,2)*rotYlm(:,:)
    rot(2*l+2:4*l+2, 1:2*l+1)=rotSpin(2,1)*rotYlm(:,:)
    rot(2*l+2:4*l+2, 2*l+2:4*l+2)=rotSpin(2,2)*rotYlm(:,:)
    !
  END SUBROUTINE

  SUBROUTINE rotate_J(rot, l, symm)
    !
    integer, intent(in)          :: l
    complex(dp), dimension(4*l+2, 4*l+2), intent(out)   :: rot
    TYPE(symmetry), intent(in)   :: symm
    !
    complex(dp), dimension(4*l+2, 4*l+2) :: Ylm2J, rot_tmp
    !
    call rotate_Ylms(rot, l, symm)
    CALL generate_Ylms2JBasis(Ylm2J, l)
    !
    !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    call zgemm('N', 'N', 4*l+2, 4*l+2, 4*l+2, cmplx_1, rot, 4*l+2, Ylm2J, 4*l+2, cmplx_0, rot_tmp, 4*l+2)
    call zgemm('C', 'N', 4*l+2, 4*l+2, 4*l+2, cmplx_1, Ylm2J, 4*l+2, rot_tmp, 4*l+2, cmplx_0, rot, 4*l+2)
    !rot=hermitian_phys_op(Ylm2J)*rot_tmp*Ylm2J
    !
  END SUBROUTINE

  SUBROUTINE rotate_Ylm(rot, l, symm)
    !
    ! Rotation matrix in Ylm basis
    !  exp(-i\alpha\hat{L}\cdot\hat{n})
    !
    use linalgwrap,  only : eigen
    !
    complex(dp), dimension(2*l+1, 2*l+1), intent(out)    :: rot
    integer, intent(in)           :: l
    TYPE(symmetry), intent(in)    :: symm
    !
    complex(dp), dimension(2*l+1, 2*l+1) :: Lx, Ly, Lz
    complex(dp), dimension(2*l+1, 2*l+1) :: Umat
    complex(dp), dimension(2*l+1, 2*l+1) :: rot_tmp
    real(dp), dimension(2*l+1)    :: eig
    complex(dp), dimension(2*l+1) :: lambda
    integer ii
    !
    if ((abs(symm%theta)<eps4).or.(l.eq.0)) then
      !
      rot=cmplx_0
      if (symm%inv.and.(mod(l,2).eq.1)) then
        do ii=1, 2*l+1
          rot(ii, ii)=-cmplx_1
        enddo
      else
        do ii=1, 2*l+1
          rot(ii, ii)=cmplx_1
        enddo
      endif
      !
    else
      !
      CALL generate_Lmatrix(Lx, Ly, Lz, l)
      !
      Umat=symm%axis(1)*Lx+symm%axis(2)*Ly+symm%axis(3)*Lz
      !
      call eigen(eig, Umat, 2*l+1)
      !
      lambda(:)=exp(-cmplx_i*symm%theta*eig(:))
      !
      if (symm%inv.and.(mod(l, 2).eq.1)) lambda(:)=-lambda(:)
      !
      do ii=1, 2*l+1
        rot_tmp(:, ii)=Umat(:, ii)*lambda(ii)
      enddo
      !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      call zgemm('N', 'C', 2*l+1, 2*l+1, 2*l+1, cmplx_1, rot_tmp, 2*l+1, Umat, 2*l+1, cmplx_0, rot, 2*l+1)
      !rot=Umat*lambda*hermitian_phys_op(Umat)
      !
    endif
    !
  END SUBROUTINE

  SUBROUTINE rotate_spinor(rot, symm)
    !
    ! Rotation matrix for spinors
    !  exp(-i\alpha\sigma\cdot\hat{n})
    !
    use linalgwrap,  only : eigen
    !
    complex(dp), dimension(2, 2), intent(out)    :: rot
    TYPE(symmetry), intent(in)    :: symm
    !
    complex(dp), dimension(2, 2) :: sx, sy, sz
    complex(dp), dimension(2, 2) :: rot_tmp
    complex(dp), dimension(2, 2) :: Umat
    real(dp), dimension(2)    :: eig
    complex(dp), dimension(2) :: lambda
    !
    integer ii
    !
    CALL generate_Smatrix(sx, sy, sz)
    !
    Umat=symm%axis(1)*sx+symm%axis(2)*sy+symm%axis(3)*sz
    !
    call eigen(eig, Umat, 2)
    !
    lambda(:)=exp(-cmplx_i*symm%theta*eig(:)/2.d0)
    !
    do ii=1, 2
      rot_tmp(:, ii)=Umat(:, ii)*lambda(ii)
    enddo
    !
    call zgemm('N', 'C', 2, 2, 2, cmplx_1, rot_tmp, 2, Umat, 2, cmplx_0, rot, 2)
    !rot=Umat*lambda*hermitian_phys_op(Umat)
    !
  END SUBROUTINE

  SUBROUTINE rotate_cubic(rot, l, symm)
    !
    ! Rotation matrix in cubic harmonic basis
    ! <C'|R|C>=<C'|Ylm'><Ylm'|R|Ylm><Ylm|C>
    !
    use linalgwrap,  only : eigen
    !
    complex(dp), dimension(2*l+1, 2*l+1), intent(out)    :: rot
    integer, intent(in)           :: l
    TYPE(symmetry), intent(in)    :: symm
    !
    complex(dp), dimension(2*l+1, 2*l+1) :: rot_ylm, ylm2c, rot_tmp
    !
    if (l.eq.0) then
      rot=1.d0
    else
      CALL generate_Ylm2C(ylm2c, l)
      CALL rotate_Ylm(rot_ylm, l, symm)
      !call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      call zgemm('N', 'N', 2*l+1, 2*l+1, 2*l+1, cmplx_1, rot_ylm, 2*l+1, ylm2c, 2*l+1, cmplx_0, rot_tmp, 2*l+1)
      call zgemm('C', 'N', 2*l+1, 2*l+1, 2*l+1, cmplx_1, ylm2c, 2*l+1, rot_tmp, 2*l+1, cmplx_0, rot, 2*l+1)
      !rot=inverse_phys_op(ylm2c)*rot_ylm*ylm2c
    endif
    !
  END SUBROUTINE

  SUBROUTINE init_symm_matrix(symm, r, t, avec)
    !
    use linalgwrap,   only: invmat
    !
    type(symmetry), intent(inout)         :: symm
    ! symmetry operator
    real(dp), dimension(3, 3), intent(in) :: r
    ! Rotation matrix (in cartesian unit) 
    real(dp), dimension(3), intent(in)    :: t
    ! Translation vector (in cartesian unit) 
    real(dp), dimension(3, 3), intent(in), optional :: avec
    ! Lattice vector (in cartesian unit)
    !
    real(dp), dimension(3, 3) :: bvec
    real(dp)  :: det
    real(dp), dimension(3)    :: wr, wi
    real(dp), dimension(12)   :: work
    real(dp), dimension(3, 3) :: tmp, vr, rt
    integer                   :: ii
    !
    if (present(avec)) then
      do ii=1, 3
        bvec(:, ii)=avec(ii, :)
      enddo
      !
      call invmat(bvec, 3)
      !
      CALL dgemm('N', 'N', 3, 3, 3, 1.d0, r, 3, avec, 3, 0.d0, tmp, 3)
      CALL dgemm('T', 'N', 3, 3, 3, 1.d0, bvec, 3, tmp, 3, 0.d0, symm%rot, 3)
    else
      symm%rot=r
    endif
    !
    symm%tau=t
    !
    det=r(1,1)*r(2,2)*r(3,3)+r(1,2)*r(2,3)*r(3,1)+r(1,3)*r(2,1)*r(3,2)- &
        r(1,3)*r(2,2)*r(3,1)-r(1,1)*r(2,3)*r(3,2)-r(1,2)*r(2,1)*r(3,3)
    symm%inv=(det<0.d0)
    !
    if (symm%inv) then
      rt(:,:)=-r(:,:)
    else
      rt(:,:)=r(:,:)
    endif
    !
    det=(rt(1,1)+rt(2,2)+rt(3,3)-1.d0)/2.d0
    if ((1.d0-det)>0.001) then
      symm%theta=acos(det)
      CALL dgeev('N', 'V', 3, rt, 3, wr, wi, tmp, 3, vr, 3, work, 12, ii)
      do ii=1, 3
        if (abs(wi(ii))<eps4.and.abs(wr(ii)-1.d0)<eps4) then
          symm%axis=vr(:, ii)
        endif
      enddo
    else
      symm%theta=0.d0
      symm%axis=(/0.d0, 0.d0, 1.d0/)
    endif
    !
  END SUBROUTINE

  SUBROUTINE init_symm_axis_angle_inv(symm, axis, theta, inv, t, avec)
    !
    use linalgwrap, only : invmat
    !
    type(symmetry), intent(inout)             :: symm
    real(dp), dimension(3), intent(in)        :: axis
    real(dp), intent(in)                      :: theta
    logical, intent(in)                       :: inv
    real(dp), dimension(3), intent(in)        :: t
    real(dp), dimension(3, 3), intent(in), optional :: avec
    !
    real(dp), dimension(3, 3)             :: bvec
    real(dp)    :: cth, sth
    real(dp)    :: n1, n2, n3
    !
    integer ii
    !
    symm%axis(:)=axis(:)/norm2(axis)
    symm%theta=theta
    symm%inv=inv
    !
    cth=cos(symm%theta)
    sth=sin(symm%theta)
    n1=symm%axis(1)
    n2=symm%axis(2)
    n3=symm%axis(3)
    !
    symm%rot(1,1)=cth+n1*n1*(1.d0-cth)
    symm%rot(1,2)=n1*n2*(1.d0-cth)-n3*sth
    symm%rot(1,3)=n1*n3*(1.d0-cth)+n2*sth
    symm%rot(2,1)=n1*n2*(1.d0-cth)+n3*sth
    symm%rot(2,2)=cth+n2*n2*(1.d0-cth)
    symm%rot(2,3)=n2*n3*(1.d0-cth)-n1*sth
    symm%rot(3,1)=n1*n3*(1.d0-cth)-n2*sth
    symm%rot(3,2)=n2*n3*(1.d0-cth)+n1*sth
    symm%rot(3,3)=cth+n3*n3*(1.d0-cth)
    !
    if (present(avec)) then
      do ii=1, 3
        bvec(:, ii)=avec(ii, :)
      enddo
      !
      CALL invmat(bvec, 3)
      !
    endif
    !
    symm%tau=t
    !
    if (inv) symm%rot=-symm%rot
    !
  END SUBROUTINE

  SUBROUTINE output_symm(symm, simple)
    !
    type(symmetry)  :: symm
    logical         :: simple
    !
    integer ii
    !
    write(stdout, '(1F7.2,A,3F7.4)', advance='no') symm%theta*360.d0/twopi, " degree rotation around axis: ", symm%axis
    !
    if (symm%inv) then
      write(stdout, *) "with inversion"
    else
      write(stdout, *)
    endif
    !
    if (.not. simple) then
      write(stdout, *) " Rotation Matrix:"
      do ii=1, 3
        write(stdout, '(3F8.3)') symm%rot(ii, :)
      enddo
      !
      write(stdout, '(A,3F8.3)') " + translation: ", symm%tau
    endif
    !
  END SUBROUTINE

  SUBROUTINE generate_Ylms2JBasis(rot, l)
    ! < Ylm s | J Jz > 
    !
    complex(dp), dimension(4*l+2, 4*l+2), intent(out)          :: rot
    integer                             :: l
    !
    integer       :: ii
    !
    rot=cmplx_0
    !
    do ii=1, 2*l          ! j=l-1/2
      ! ii=   1  ,    2  , ... 2l-1 , 2l
      ! mj=-l+1/2, -l+3/2, ... l-3/2, l-1/2  
      ! == For <lm, up|, mj=m+1/2
      ! m=ii-(l+1)
      ! jj=m+l+1          ! <lm+| : <lm,up|
      rot(ii, ii)=sqrt((2*l+1-ii)/(2.d0*l+1.d0))
      ! == For <lm, dn|, mj=m-1/2
      ! m=ii-l
      ! jj=m+3*l+2        ! <lm-| : <lm,dn|
      rot(ii+2*l+2, ii)=-sqrt(ii/(2.d0*l+1.d0))
    enddo
    !
    do ii=2*l+1, 4*l+2    ! j=l+1/2
      ! ii=  2l+1,   2l+2, ...  4l+1,  4l+2
      ! mj=-l-1/2, -l+1/2, ... l-1/2, l+1/2  
      ! == For <lm, up|, mj=m+1/2
      ! m=ii-3*l-2
      ! jj=m+l+1          ! <lm+| : <lm,up|
      if (ii>2*l+1) rot(ii-2*l-1, ii)=sqrt((ii-2.d0*l-1)/(2.d0*l+1.d0))
      ! == For <lm, dn|, mj=m-1/2
      ! m=ii-3*l-1
      ! jj=m+3*l+2        ! <lm-| : <lm,dn|
      if (ii<4*l+2) rot(ii+1, ii)=sqrt((4*l+2.d0-ii)/(2.d0*l+1.d0))
    enddo
    !
  END SUBROUTINE
    !
END MODULE
