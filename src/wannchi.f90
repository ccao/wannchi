include 'lapack.f90'

PROGRAM wannchi
  !
  use mpi
  use constants
  use wanndata
  use estate_mesh
  !
  implicit none
  !
  integer iproc, nproc
  integer first_q, last_q
  integer nqmesh
  integer nqx, nqy, nqz
  integer iqx, iqy, iqz, iq
  real(dp) qvec(1:3)
  complex(dp) chi
  character(len=80) seed
  !
  seed="wannier90"
  !
  CALL read_ham(seed)
  !
  open(unit=fin, file="input")
  read(fin, *) ef
  read(fin, *) nkx, nky, nkz
  close(unit=fin)
  !
  nkmesh=nkx*nky*nkz
  !
  nqx=nkx; nqy=nky; nqz=1
  nqmesh=nqx*nqy*nqz
  !
  allocate(eig(1:nkmesh, 1:norb))
  allocate(egv(1:nkmesh, 1:norb, 1:norb))
  !
  CALL interpolate_states
  !
!  do iq=1, nqmesh
  do iqx=1, nqx
    do iqy=1, nqy
    qvec(1)=(iqx-1.d0)/nqx
    qvec(2)=(iqy-1.d0)/nqy
    qvec(3)=0.d0
    CALL compute_chi_diag(chi, qvec)
    !
    write(*,'(3F12.8,2F22.12)') qvec, chi
    enddo
    write(*,*) ' '
  enddo
!  enddo
  !
  deallocate(eig, egv)
  !
  CALL finalize_wann
  !
END PROGRAM

SUBROUTINE interpolate_states
  !
  use lapack95
  use constants
  use wanndata
  use estate_mesh
  !
  implicit none
  !
  real(dp) kvec(1:3)
  real(dp) rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:,:)
  real(dp), allocatable :: e(:)
  !
  integer ix, iy, iz, ir, ik, info
  !
  allocate(work(1:norb, 1:norb))
  allocate(e(1:norb))
  !
  write(stdout, *) "Starting interpolation of the original states to"
  write(stdout, *) nkx, "x", nky, "x", nkz, " K-mesh"
  do ix=1, nkx
    kvec(1)=(ix-1.d0)/nkx
    do iy=1, nky
      kvec(2)=(iy-1.d0)/nky
      do iz=1, nkz
        kvec(3)=(iz-1.d0)/nkz
        ik=(ix-1)*nky*nkz+(iy-1)*nky+iz
        !
        work(:,:)=cmplx_0
        do ir=1, nrpt
          rdotk=SUM(kvec(:)*rvec(ir, :))
          fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
          work(:,:)=work(:,:)+ham(:,:,ir)*fact
        enddo ! ir
        call heev(work, e, 'V', 'U', info)
        !
        egv(ik, :, :)=work(:,:)
        eig(ik, :)=e(:)
      enddo ! iz
    enddo ! iy
  enddo ! ix
  write(stdout, *) "Done..."
  !
  deallocate(work)
  !
END SUBROUTINE

FUNCTION norm(a, b)
  !
  implicit none
  !
  integer norm, a, b
  !
  norm=modulo(a, b)
  if (norm.eq.0) norm=b
  !
  return
END FUNCTION

SUBROUTINE compute_chi_diag(chi, qv)
  !
  use constants
  use wanndata
  use estate_mesh
  !
  implicit none
  !
  complex(dp) chi
  !
  real(dp) qv(1:3)
  !
  integer, external :: norm
  !
  integer ix, iy, iz, jx, jy, jz, dx, dy, dz, ik, jk
  integer io, jo, mu, nu
  complex(dp) fact
  real(dp) eta
  !
  eta=1.0d-3
  !
  dx=qv(1)*nkx; dy=qv(2)*nky; dz=qv(3)*nkz
  chi=cmplx_0
  !
  do ix=1, nkx
    do iy=1, nky
      do iz=1, nkz
        jx=norm(ix+dx, nkx); jy=norm(iy+dy, nky); jz=norm(iz+dz, nkz)
        ik=(ix-1)*nky*nkz+(iy-1)*nkz+iz
        jk=(jx-1)*nky*nkz+(jy-1)*nkz+jz
        !
        do io=1, norb
          do jo=1, norb
            if ( (eig(ik, io)>ef).and.(eig(jk, jo)<ef) ) then
              fact=SUM(egv(ik, :, io)*conjg(egv(jk, :, jo)))
              chi=chi+conjg(fact)*fact/(eig(jk, jo)-eig(ik, io)+eta*cmplx_i)
            else
              if ( (eig(ik, io)<ef).and.(eig(jk, jo)>ef) ) then
                fact=SUM(egv(ik, :, io)*conjg(egv(jk, :, jo)))
                chi=chi-conjg(fact)*fact/(eig(jk, jo)-eig(ik, io)-eta*cmplx_i)
              endif
            endif
          enddo
        enddo
        !
      enddo
    enddo
  enddo
  !
  chi=-chi/(nkx*nky*nkz)
  !
END SUBROUTINE
