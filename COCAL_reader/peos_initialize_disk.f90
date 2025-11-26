subroutine peos_initialize_disk(dir_path)
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter_disk   !abc_disk,abi_disk,rhoi_disk,qi_disk,hi_disk,nphase_disk,rhoini_cgs_disk,emdini_gcm1_disk
  implicit none
  character*400, intent(in) :: dir_path
  real(8) :: rho_0, pre_0, facrho, facpre, fac2
  integer :: ii, iphase
!
  open(852,file=trim(dir_path)//'/'//'peos_parameter_NStorus.dat',status='old')
  read(852,'(8x,1i5,es13.5)') nphase_disk, rhoini_cgs_disk
  read(852,'(2es13.5)') rho_0, pre_0
  do ii = nphase_disk, 0, -1
    read(852,'(2es13.5)') rhocgs_disk(ii), abi_disk(ii)
  end do
  close(852)
!
! --  cgs to g = c = msol = 1 unit.
! --  assume pre = pre_0 dyn/cm^2 at rho = rho_0 gr/cm^3.
! --  typically pre_0 = 1.0d+37 dyn/cm^2 
! --  and       rho_0 = 1.0d+16  gr/cm^3.
! --  rescale interface values
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8
!
  do ii = 0, nphase_disk
    rhoi_disk(ii) = facrho*rhocgs_disk(ii)
  end do
!
  call peos_lookup(rho_0,rhocgs_disk,nphase_disk,iphase)
!    
  abc_disk(iphase) = pre_0/rho_0**abi_disk(iphase)
  abc_disk(iphase) = facpre/facrho**abi_disk(iphase)*abc_disk(iphase)
  abccgs_disk(iphase) = pre_0/(rho_0**abi_disk(iphase))
!
  if (iphase.gt.0) then
    do ii = iphase-1, 0, -1
      abc_disk(   ii) = rhoi_disk(  ii)**(abi_disk(ii+1)-abi_disk(ii))*abc_disk(   ii+1)
      abccgs_disk(ii) = rhocgs_disk(ii)**(abi_disk(ii+1)-abi_disk(ii))*abccgs_disk(ii+1)
    end do
  end if
  if (iphase.lt.nphase_disk) then
    do ii = iphase+1, nphase_disk
      abc_disk(   ii) = rhoi_disk(  ii-1)**(abi_disk(ii-1)-abi_disk(ii))*abc_disk(   ii-1)
      abccgs_disk(ii) = rhocgs_disk(ii-1)**(abi_disk(ii-1)-abi_disk(ii))*abccgs_disk(ii-1)
    end do
  end if
!
  do ii = 0, nphase_disk
    qi_disk(ii) = abc_disk(ii)*rhoi_disk(ii)**(abi_disk(ii)-1.0d0)
  end do
!
  hi_disk(0) = 1.0d0
  do ii = 1, nphase_disk
    fac2 = abi_disk(ii)/(abi_disk(ii) - 1.0d0)
    hi_disk(ii) = hi_disk(ii-1) + fac2*(qi_disk(ii) - qi_disk(ii-1))
  end do
!
  rhoini_gcm1_disk = facrho*rhoini_cgs_disk
  call peos_lookup(rhoini_gcm1_disk,rhoi_disk,nphase_disk,iphase)
  emdini_gcm1_disk = abc_disk(iphase)*rhoini_gcm1_disk**(abi_disk(iphase)-1.0d0)
!
end subroutine peos_initialize_disk