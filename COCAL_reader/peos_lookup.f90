subroutine peos_lookup(qp,qpar,nphase_in,iphase)
!
  use phys_constant         !nnpeos
  implicit none
!
  real(8), intent(in)  :: qp, qpar(0:nnpeos)
  integer, intent(in)  :: nphase_in
  real(8)              :: det
  integer, intent(out) :: iphase
  integer              :: ii
!
! --  Monotonically increasing qpar is assumed.
!
  iphase = 1
  do ii = 1, nphase_in
    det = (qp-qpar(ii))*(qp-qpar(ii-1))
    if (det <= 0.0d0) then
      iphase = ii
      exit
    end if
  end do
!
end subroutine peos_lookup