module def_peos_parameter_NS
  use phys_constant         !nnpeos
  implicit none
  real(8) :: abc_NS(0:nnpeos), abi_NS(0:nnpeos), rhoi_NS(0:nnpeos), &
  &          qi_NS(0:nnpeos), hi_NS(0:nnpeos)
  real(8) :: rhocgs_NS(0:nnpeos), abccgs_NS(0:nnpeos)
  real(8) :: rhoini_cgs_NS, rhoini_gcm1_NS, emdini_gcm1_NS
  real(8) :: sgma_NS, constqc_NS, cbar_NS
  integer :: nphase_NS
end module def_peos_parameter_NS