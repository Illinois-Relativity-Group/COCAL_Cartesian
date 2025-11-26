subroutine IO_input_CF_grav_export(filenm,psi,alph,bvxd,bvyd,bvzd)
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer :: psi(:,:,:), alph(:,:,:), bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:)
  character(len=*) :: filenm
!
  write(6,*) "Reading psi, alpha, beta_i..."
! --- Metric potentials.
  open(13,file=trim(filenm),status='old')
  read(13,'(5i5)') nrtmp, nttmp, nptmp
  do ip = 1, nptmp
    do it = 1, nttmp
      do ir = 1, nrtmp
        !print *, 'Reading data for ir =', ir, ', it =', it, ', ip =', ip
        read(13,'(1p,5e20.12)')  psi(ir,it,ip), &
    &                           alph(ir,it,ip), &
    &                           bvxd(ir,it,ip), &
    &                           bvyd(ir,it,ip), &
    &                           bvzd(ir,it,ip)
        !print *, 'Successfully read:', psi(ir,it,ip)
      end do
    end do
  end do
  close(13)
!
end subroutine IO_input_CF_grav_export
