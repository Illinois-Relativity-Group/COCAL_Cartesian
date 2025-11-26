subroutine peos_initialize_NS(dir_path)
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter_NS   !abc_NS,abi_NS,rhoi_NS,qi_NS,hi_NS,nphase_NS,rhoini_cgs_NS,emdini_gcm1_NS
  implicit none
  character*400, intent(in) :: dir_path
  real(8) :: rho_0, pre_0, facrho, facpre, fac2
  integer :: ii, iphase
!
  open(851,file=trim(dir_path)//'/'//'peos_parameter.dat',status='old')
  read(851,'(8x,1i5,es13.5)') nphase_NS, rhoini_cgs_NS
  read(851,'(2es13.5)') rho_0, pre_0
  do ii = nphase_NS, 0, -1
    read(851,'(2es13.5)') rhocgs_NS(ii), abi_NS(ii)
  end do
  close(851)
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
  do ii = 0, nphase_NS
    rhoi_NS(ii) = facrho*rhocgs_NS(ii)
  end do
!
  call peos_lookup(rho_0,rhocgs_NS,nphase_NS,iphase)
!    
  abc_NS(iphase) = pre_0/rho_0**abi_NS(iphase)
  abc_NS(iphase) = facpre/facrho**abi_NS(iphase)*abc_NS(iphase)
  abccgs_NS(iphase) = pre_0/(rho_0**abi_NS(iphase))
!
  if (iphase.gt.0) then
    do ii = iphase-1, 0, -1
      abc_NS(   ii) = rhoi_NS(  ii)**(abi_NS(ii+1)-abi_NS(ii))*abc_NS(   ii+1)
      abccgs_NS(ii) = rhocgs_NS(ii)**(abi_NS(ii+1)-abi_NS(ii))*abccgs_NS(ii+1)
    end do
  end if
  if (iphase.lt.nphase_NS) then
    do ii = iphase+1, nphase_NS
      abc_NS(   ii) = rhoi_NS(  ii-1)**(abi_NS(ii-1)-abi_NS(ii))*abc_NS(   ii-1)
      abccgs_NS(ii) = rhocgs_NS(ii-1)**(abi_NS(ii-1)-abi_NS(ii))*abccgs_NS(ii-1)
    end do
  end if
!
  do ii = 0, nphase_NS
    qi_NS(ii) = abc_NS(ii)*rhoi_NS(ii)**(abi_NS(ii)-1.0d0)
  end do
!
  hi_NS(0) = 1.0d0
  do ii = 1, nphase_NS
    fac2 = abi_NS(ii)/(abi_NS(ii) - 1.0d0)
    hi_NS(ii) = hi_NS(ii-1) + fac2*(qi_NS(ii) - qi_NS(ii-1))
  end do
!
  rhoini_gcm1_NS = facrho*rhoini_cgs_NS
  call peos_lookup(rhoini_gcm1_NS,rhoi_NS,nphase_NS,iphase)
  emdini_gcm1_NS = abc_NS(iphase)*rhoini_gcm1_NS**(abi_NS(iphase)-1.0d0)
!
end subroutine peos_initialize_NS