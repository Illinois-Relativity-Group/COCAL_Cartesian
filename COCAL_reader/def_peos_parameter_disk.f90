module def_peos_parameter_disk
  use phys_constant         !nnpeos
  implicit none
  real(8) :: abc_disk(0:nnpeos), abi_disk(0:nnpeos), rhoi_disk(0:nnpeos), &
  &          qi_disk(0:nnpeos), hi_disk(0:nnpeos)
  real(8) :: rhocgs_disk(0:nnpeos), abccgs_disk(0:nnpeos)
  real(8) :: rhoini_cgs_disk, rhoini_gcm1_disk, emdini_gcm1_disk
  real(8) :: sgma_disk, constqc_disk, cbar_disk
  integer :: nphase_disk
end module def_peos_parameter_disk