subroutine peos_q2hprho_disk(q,h,pre,rho,ened)
!
  use def_peos_parameter_disk        !abc_disk,abi_disk,rhoi_disk,qi_disk,hi_disk,nphase_disk
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: h, pre, rho
  real(8)                :: hin, qin, abin, abct, fac1, fac2, fack, small, ened
  integer                :: iphase
!
  call peos_lookup(q, qi_disk, nphase_disk, iphase)
  hin  = hi_disk(iphase)
  qin  = qi_disk(iphase)
  abin = abi_disk(iphase)
  abct = abc_disk(iphase)
!
  fac1 = 1.0d0/(abin - 1.0d0)
  fac2 = abin/(abin - 1.0d0)
  fack = abct**(-fac1)
!
  small = 1.0d-60
  if (q <= small) q = small
  h = hin + fac2*(q - qin)
  if (h <= 1.0d0) h = 1.0d0
  pre = fack*q**fac2
  rho = fack*q**fac1
  ened = rho*h - pre
!
end subroutine peos_q2hprho_disk