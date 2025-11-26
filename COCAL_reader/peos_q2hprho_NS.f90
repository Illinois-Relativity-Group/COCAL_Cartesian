subroutine peos_q2hprho_NS(q,h,pre,rho,ened)
!
  use def_peos_parameter_NS        !abc_NS,abi_NS,rhoi_NS,qi_NS,hi_NS,nphase_NS
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: h, pre, rho
  real(8)                :: hin, qin, abin, abct, fac1, fac2, fack, small, ened
  integer                :: iphase
!
  call peos_lookup(q, qi_NS, nphase_NS, iphase)
  hin  = hi_NS(iphase)
  qin  = qi_NS(iphase)
  abin = abi_NS(iphase)
  abct = abc_NS(iphase)
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
end subroutine peos_q2hprho_NS