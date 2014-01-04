INTEGER FUNCTION rpt_translated(irpt, drpt)
  !
  use constants
  use wanndata, only: rvec, nrpt
  !
  IMPLICIT NONE
  !
  integer irpt, drpt(1:3)
  integer new_rvec(1:3), ii
  logical found
  !
  found=.false.
  new_rvec(:)=rvec(:,irpt)+drpt(:)
  do ii=1, nrpt
    if ((new_rvec(1).eq.rvec(1,ii)).and.(new_rvec(2).eq.rvec(2,ii)).and.(new_rvec(3).eq.rvec(3,ii))) then
      found=.true.
      rpt_translated=ii
      return
    endif
  enddo
  if (.not.found) rpt_translated=-1
  !
END FUNCTION

SUBROUTINE mix_ham
  !
  use constants
  use wanndata
  use orbitals
  use dopedata
  !
  IMPLICIT NONE
  !
  integer iatm, jatm, iiatm, jjatm
  integer iorb, jorb, iiorb, jjorb
  integer irpt, iirpt, ii, rpt_translated
  real(dp) trans_vec(1:3)    ! : denotes the relative primitive cell shift of dopant to reference
  integer drpt(1:3)
  !
  doped_ham(:,:,:)=ham(0,:,:,:)
  do ii=1, ndopant
    write(stdout, *) " !!!: Dopant ", ii
    write(stdout, *) " Coordinates: "
    write(stdout, *)  tau(dopant_list(ii), :)
    trans_vec(:)=tau(dopant_list(ii), :)-tau(ref_dopant, :)
    CALL translate_crystal(trans_vec)
    do iatm=1, nat
      if (sp_norb(at_sp(iatm)).eq.0) cycle
      do jatm=1, nat
        if (sp_norb(at_sp(jatm)).eq.0) cycle
        iiatm=new_to_old(iatm)
        jjatm=new_to_old(jatm)
        drpt(:)=new_to_old_dR(jatm, :)-new_to_old_dR(iatm, :)
        do irpt=1, nrpt
          iirpt=rpt_translated(irpt, drpt)
          do iorb=1, sp_norb(at_sp(iatm))
            iiorb=glob_orb(iiatm, iorb)
            do jorb=1, sp_norb(at_sp(jatm))
              jjorb=glob_orb(jjatm, jorb)
              if (iirpt.ne.-1) then
                dham(glob_orb(iatm, iorb), glob_orb(jatm, jorb), irpt)=ham(1,iiorb, jjorb, iirpt)
              else
                dham(glob_orb(iatm, iorb), glob_orb(jatm, jorb), irpt)=cmplx_0
              endif
            enddo
          enddo ! iorb
        enddo ! irpt
      enddo
    enddo
    doped_ham(:,:,:)=doped_ham(:,:,:)+dham(:,:,:)
  enddo
  !
END SUBROUTINE
