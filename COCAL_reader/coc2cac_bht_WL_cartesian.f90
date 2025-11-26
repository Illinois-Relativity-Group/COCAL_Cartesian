!_____________________________________________________________________________
!
!    CACTUS READER OF COCAL BHT ID IN WAVELESS FORMALISM 
!_____________________________________________________________________________
!
include './phys_constant.f90'
include './def_matter_parameter.f90'
include './def_quantities.f90'
include './def_bh_parameter.f90'
include './def_bht_parameter.f90'
include './def_peos_parameter.f90'
include './make_array_2d.f90'
include './make_array_3d.f90'
include './grid_parameter.f90'
include './interface_modules_cartesian.f90' 
include './coordinate_grav_r.f90' 
include './coordinate_grav_phi.f90' 
include './coordinate_grav_theta.f90' 
include './coordinate_grav_extended.f90' 
include './trigonometry_grav_theta.f90' 
include './trigonometry_grav_phi.f90' 
include './calc_bht_excision_radius.f90' 
include './read_bht_parameter.f90' 
include './interface_IO_input_CF_grav_export.f90' 
include './interface_IO_input_WL_grav_export_hij.f90' 
include './interface_IO_input_grav_export_Kij.f90' 
include './interface_IO_input_CF_star_export.f90' 
include './interface_IO_input_matter_BHT_export.f90' 
include './interface_invhij_WL_export.f90' 
include './interface_index_vec_down2up_export.f90' 
include './interface_interpo_gr2fl_metric_CF_export.f90' 
include './interface_IO_input_grav_export_Ai.f90' 
include './interface_IO_input_grav_export_Faraday.f90' 
include './interface_IO_input_star4ve_export.f90' 
include './interface_interpo_gr2fl_export.f90' 
include './interface_interpo_lag4th_2Dsurf.f90' 
include './IO_input_CF_grav_export.f90'
include './IO_input_WL_grav_export_hij.f90'
include './IO_input_grav_export_Kij.f90'
include './IO_input_CF_star_export.f90'
include './IO_input_matter_BHT_export.f90' 
include './invhij_WL_export.f90'
include './index_vec_down2up_export.f90'
include './interpo_gr2fl_metric_CF_export.f90'
include './IO_input_grav_export_Ai.f90'
include './IO_input_grav_export_Faraday.f90'
include './IO_input_star4ve_export.f90'
include './interpo_gr2fl_export.f90'
include './interpo_gr2cgr_4th.f90'
include './interpo_fl2cgr_4th_export.f90'
include './interpo_lag4th_2Dsurf.f90' 
include './lagint_4th.f90' 
include './peos_initialize.f90' 
include './peos_q2hprho.f90' 
include './peos_lookup.f90' 
!
!___________________________________________________

PROGRAM coc2cac
  use phys_constant
  use grid_parameter
  use interface_modules_cartesian
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use interface_IO_input_CF_grav_export
  use interface_IO_input_WL_grav_export_hij
  use interface_IO_input_grav_export_Kij
  use interface_IO_input_matter_BHT_export
  use interface_invhij_WL_export
  use interface_index_vec_down2up_export
  implicit none
  integer :: iAB, ix, iy, iz, res, i, file_unit
  character(30) :: char1, filename, temp_str
  character*400 :: dir_path
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc, rcoc
  real(8) :: emdgca, omegca
  real(8) :: psica,  alphca, psi4ca
  real(8) :: bvxdca, bvydca, bvzdca, bvxuca, bvyuca, bvzuca
  real(8) :: hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca
  real(8) :: hxxuca, hxyuca, hxzuca, hyyuca, hyzuca, hzzuca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: kxxca, kxyca, kxzca, kyyca, kyzca, kzzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor
  real(8) :: gxx1, gxy1, gxz1, gyy1, gyz1, gzz1, kxx1, kxy1, kxz1, kyy1, kyz1, kzz1
  real(8) :: ome, ber, radi, rexc

  real(long) :: x_min, x_max, y_min, y_max, z_min, z_max
  real(long) :: dx, dy, dz
!
  real(8), pointer :: emdg(:,:,:) , omeg(:,:,:)
  real(8), pointer :: psi(:,:,:)  , alph(:,:,:)
  real(8), pointer :: bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:) , bvxu(:,:,:) , bvyu(:,:,:), bvzu(:,:,:)
  real(8), pointer :: hxxd(:,:,:) , hxyd(:,:,:) , hxzd(:,:,:) , hyyd(:,:,:) , hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:) , hxyu(:,:,:) , hxzu(:,:,:) , hyyu(:,:,:) , hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: kxx(:,:,:),   kxy(:,:,:)  , kxz(:,:,:)  , kyy(:,:,:)  , kyz(:,:,:) , kzz(:,:,:)

  real(long), pointer :: x_cart(:), y_cart(:), z_cart(:)

  !real(long), pointer :: rho_arr(:,:,:)
!
    call get_command_argument(1, temp_str)
    read(temp_str,*) x_min

    call get_command_argument(2, temp_str)
    read(temp_str,*) x_max

    call get_command_argument(3, temp_str)
    read(temp_str,*) y_min

    call get_command_argument(4, temp_str)
    read(temp_str,*) y_max

    call get_command_argument(5, temp_str)
    read(temp_str,*) z_min

    call get_command_argument(6, temp_str)
    read(temp_str,*) z_max

    call get_command_argument(7, temp_str)
    read(temp_str,*) res

    allocate(x_cart(res))
    allocate(y_cart(res))
    allocate(z_cart(res))

    dx = (x_max - x_min) / (res - 1)
    dy = (y_max - y_min) / (res - 1)
    dz = (z_max - z_min) / (res - 1)

    do i = 1, res
      x_cart(i) = x_min + (i - 1) * dx
      y_cart(i) = y_min + (i - 1) * dy
      z_cart(i) = z_min + (i - 1) * dz
    end do


  gxx1=0.0d0; gxy1=0.0d0; gxz1=0.0d0; gyy1=0.0d0; gyz1=0.0d0; gzz1=0.0d0
  kxx1=0.0d0; kxy1=0.0d0; kxz1=0.0d0; kyy1=0.0d0; kyz1=0.0d0; kzz1=0.0d0
  kxxca=0.0d0; kxyca=0.0d0; kxzca=0.0d0; kyyca=0.0d0; kyzca=0.0d0; kzzca=0.0d0

!  dir_path=dir_path(1:path_len)
!  write(*,*) "In coc2pri_bht_WL: dir_path=",dir_path

  write(*,*) "Reading initial data from..."
  open(70,file='PATH2ID.txt',status='old')
  read(70,'(i1)') iAB
  read(70,'(a400)') dir_path
  close(70)
  write(*,*) dir_path 

! -- Read parameters
  call read_parameter_cactus(dir_path)
  call read_bht_parameter_cactus(dir_path)
  call calc_bht_excision_radius
  call peos_initialize_cactus(dir_path)
  call grid_r_bht('eBH')
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call grid_extended
!
  rexc = rg(0)
  write(6,'(a14,1p,1e23.15)') "Excision at r=", rexc
!    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
!  rr3 = 0.7d0*(rgout - rgmid)
!  dis_cm = dis

  allocate ( emdg(0:nrg,0:ntg,0:npg))
  allocate ( omeg(0:nrg,0:ntg,0:npg))
  allocate (  psi(0:nrg,0:ntg,0:npg))
  allocate ( alph(0:nrg,0:ntg,0:npg))
  allocate ( bvxd(0:nrg,0:ntg,0:npg))
  allocate ( bvyd(0:nrg,0:ntg,0:npg))
  allocate ( bvzd(0:nrg,0:ntg,0:npg))
  allocate ( bvxu(0:nrg,0:ntg,0:npg))
  allocate ( bvyu(0:nrg,0:ntg,0:npg))
  allocate ( bvzu(0:nrg,0:ntg,0:npg))
  allocate ( hxxd(0:nrg,0:ntg,0:npg))
  allocate ( hxyd(0:nrg,0:ntg,0:npg))
  allocate ( hxzd(0:nrg,0:ntg,0:npg))
  allocate ( hyyd(0:nrg,0:ntg,0:npg))
  allocate ( hyzd(0:nrg,0:ntg,0:npg))
  allocate ( hzzd(0:nrg,0:ntg,0:npg))
  allocate ( hxxu(0:nrg,0:ntg,0:npg))
  allocate ( hxyu(0:nrg,0:ntg,0:npg))
  allocate ( hxzu(0:nrg,0:ntg,0:npg))
  allocate ( hyyu(0:nrg,0:ntg,0:npg))
  allocate ( hyzu(0:nrg,0:ntg,0:npg))
  allocate ( hzzu(0:nrg,0:ntg,0:npg))
  allocate (  kxx(0:nrg,0:ntg,0:npg))
  allocate (  kxy(0:nrg,0:ntg,0:npg))
  allocate (  kxz(0:nrg,0:ntg,0:npg))
  allocate (  kyy(0:nrg,0:ntg,0:npg))
  allocate (  kyz(0:nrg,0:ntg,0:npg))
  allocate (  kzz(0:nrg,0:ntg,0:npg))

  !allocate(rho_arr(0:res,0:res,0:res))

  emdg=0.0d0; omeg=0.0d0
  psi=0.0d0;  alph=0.0d0
  bvxd=0.0d0; bvyd=0.0d0;  bvzd=0.0d0;  bvxu=0.0d0;  bvyu=0.0d0;  bvzu=0.0d0
  kxx=0.0d0;  kxy =0.0d0;  kxz =0.0d0;   kyy=0.0d0;   kyz=0.0d0;   kzz=0.0d0
  hxxd=0.0d0; hxyd=0.0d0;  hxzd=0.0d0;  hyyd=0.0d0;  hyzd=0.0d0;  hzzd=0.0d0;
  hxxu=0.0d0; hxyu=0.0d0;  hxzu=0.0d0;  hyyu=0.0d0;  hyzu=0.0d0;  hzzu=0.0d0;


  call IO_input_CF_grav_export(trim(dir_path)//'/'//'rnsgra_3D.las',psi,alph,bvxd,bvyd,bvzd)

  call IO_input_WL_grav_export_hij(trim(dir_path)//'/'//'rnsgra_hij_3D.las',hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)

  call IO_input_grav_export_Kij(trim(dir_path)//'/'//'rnsgra_Kij_3D.las',kxx,kxy,kxz,kyy,kyz,kzz)

  call IO_input_matter_BHT_export(trim(dir_path)//"/rnsflu_3D.las",emdg,omeg,ome,ber,radi)

  call invhij_WL_export(hxxd,hxyd,hxzd,hyyd,hyzd,hzzd,hxxu,hxyu,hxzu,hyyu,hyzu,hzzu)

  call index_vec_down2up_export(hxxu,hxyu,hxzu,hyyu,hyzu,hzzu,bvxu,bvyu,bvzu,bvxd,bvyd,bvzd)


  write(6,'(2e20.12)') emdg(0,0,0), omeg(0,0,0)
  write(6,'(3e20.12,1I4.2)') ome, ber, radi

  filename = "data.txt"

  open(unit=file_unit, file=filename, status='replace', action='write')
  write(file_unit, '(A23, A23, A23, A23, A23, A23, A23, A23, & 
                  A23, A23, A23, A23, A23, A23, A23, A23, & 
                  A23, A23, A23, A23, A23, A23, A23, A23, A23, A23, A23, A23, A23)') &
  'x', 'y', 'z', 'rho', 'psi', 'alph', 'bvxd', 'bvyd', 'bvzd', 'bvxu', &
  'bvyu', 'bvzu', 'hxxd', 'hxyd', 'hxzd', 'hyyd', 'hyzd', 'hzzd', 'kxx', &
  'kxy', 'kxz', 'kyy', 'kyz', 'kzz', 'emdg', 'omeg', 'h', 'ene', 'eps'


!  
do ix = 1, res
  do iy = 1, res
    do iz = 1, res
      !print *, "Point (", x_cart(ix), ", ", y_cart(iy), ", ", z_cart(iz), ")"

      !write(6,'(a56)', ADVANCE = "NO") "Give cartesian coordinates (x,y,z) separated by a space:"
      !read(5,*) xcac,ycac,zcac
      xcac = x_cart(ix)
      ycac = y_cart(iy)
      zcac = z_cart(iz)
      !write(6,'(a23,3e20.12)') "Point given wrt CACTUS:", xcac,ycac,zcac
      !write(6,'(a20,1e20.12)') "Cocal radius scale :", radi
      xcoc = xcac/(radi)
      ycoc = ycac/(radi)
      zcoc = zcac/(radi)
      !write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc


      rcoc = dsqrt(xcoc*xcoc + ycoc*ycoc + zcoc*zcoc)

    !  if (rcoc < 1.1d0*rexc)  cycle
      if (rcoc < rexc) then
        write(6,*) 'Point is inside excised radius. Exiting...'
        stop  
      end if

      call interpo_gr2cgr_4th(psi , psica , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(alph, alphca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvxd, bvxdca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvyd, bvydca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvzd, bvzdca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvxu, bvxuca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvyu, bvyuca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(bvzu, bvzuca, xcoc, ycoc, zcoc)

      call interpo_gr2cgr_4th(hxxd, hxxdca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(hxyd, hxydca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(hxzd, hxzdca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(hyyd, hyydca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(hyzd, hyzdca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(hzzd, hzzdca, xcoc, ycoc, zcoc)

      call interpo_gr2cgr_4th(kxx , kxxca , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(kxy , kxyca , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(kxz , kxzca , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(kyy , kyyca , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(kyz , kyzca , xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(kzz , kzzca , xcoc, ycoc, zcoc)

      call interpo_gr2cgr_4th(emdg, emdgca, xcoc, ycoc, zcoc)
      call interpo_gr2cgr_4th(omeg, omegca, xcoc, ycoc, zcoc)

      bxcor = bvxuca + omegca*(-ycoc)
      bycor = bvyuca + omegca*(xcoc)
      bzcor = bvzuca
      psi4ca = psica**4

      if (dabs(emdgca) > 1.0d-14) then
        vxu = bxcor/alphca 
        vyu = bycor/alphca
        vzu = bzcor/alphca
      else
        emdgca=0.0d0
        vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
      end if

      ! Setting 3-metric and assigning Kij to gbij_t
      gxx1 = psi4ca*(1.0d0+hxxdca)
      gxy1 = psi4ca*(      hxydca)
      gxz1 = psi4ca*(      hxzdca)
      gyy1 = psi4ca*(1.0d0+hyydca)
      gyz1 = psi4ca*(      hyzdca)
      gzz1 = psi4ca*(1.0d0+hzzdca)

      kxx1 = psi4ca*kxxca/(radi)
      kxy1 = psi4ca*kxyca/(radi)
      kxz1 = psi4ca*kxzca/(radi)
      kyy1 = psi4ca*kyyca/(radi)
      kyz1 = psi4ca*kyzca/(radi)
      kzz1 = psi4ca*kzzca/(radi)

    !  beta_x = psi4ca*bvxdca     ! b_x lower index
    !  beta_y = psi4ca*bvydca
    !  beta_z = psi4ca*bvzdca
      
      call peos_q2hprho(emdgca, hca, preca, rhoca, eneca)


      epsca = eneca/rhoca - 1.0d0

      write(file_unit, '(E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, & 
    & E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, & 
    & E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, E23.15, & 
    & E23.15, E23.15, E23.15, E23.15, E23.15)') &
      xcac, ycac, zcac, rhoca, psica, alphca, bvxdca, bvydca, bvzdca, bvxuca, &
      bvyuca, bvzuca, hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca, kxxca, &
      kxyca, kxzca, kyyca, kyzca, kzzca, emdgca, omegca, hca, eneca, epsca


      !write(6,'(a6,e20.12)') "psi  =", psica
      !write(6,'(a6,e20.12)') "alph =", alphca
      !write(6,'(a6,e20.12)') "bvxd =", bvxdca
      !write(6,'(a6,e20.12)') "bvyd =", bvydca
      !write(6,'(a6,e20.12)') "bvzd =", bvzdca
      !write(6,'(a6,e20.12)') "Radi =", radi
      !write(6,'(a6,e20.12)') "Omeg =", ome/radi
      !write(6,'(a6,e20.12)') "emdg =", emdgca
      !write(6,'(a6,e20.12)') "h    =", hca
      !write(6,'(a6,e20.12)') "pre  =", preca
      !write(6,'(a6,e20.12)') "rho  =", rhoca
      !write(6,'(a6,e20.12)') "ene  =", eneca
      !write(6,'(a6,e20.12)') "eps  =", epsca
      !
      !write(6,'(a18)') "gij at gridpoints:"
      !write(6,'(3e20.12)') gxx1, gxy1, gxz1
      !write(6,'(3e20.12)') gxy1, gyy1, gyz1
      !write(6,'(3e20.12)') gxz1, gyz1, gzz1

      !write(6,'(a18)') "Kij at gridpoints:"
      !write(6,'(3e20.12)') kxx1, kxy1, kxz1
      !write(6,'(3e20.12)') kxy1, kyy1, kyz1
      !write(6,'(3e20.12)') kxz1, kyz1, kzz1

      !write(6,'(a13)') "v^i eulerian:"
      !write(6,'(a6,e20.12)') "vxu  =", vxu
      !write(6,'(a6,e20.12)') "vyu  =", vyu
      !write(6,'(a6,e20.12)') "vzu  =", vzu
    end do
  end do
end do
           
  write(6,'(a16)') "Deallocating...."
  deallocate( emdg);  deallocate( omeg);
  deallocate(  psi);  deallocate( alph);  deallocate( bvxd);  deallocate( bvyd);
  deallocate( bvzd);  deallocate( bvxu);  deallocate( bvyu);  deallocate( bvzu);
  deallocate( hxxd);  deallocate( hxyd);  deallocate( hxzd);  deallocate( hyyd);
  deallocate( hyzd);  deallocate( hzzd);  deallocate( hxxu);  deallocate( hxyu);
  deallocate( hxzu);  deallocate( hyyu);  deallocate( hyzu);  deallocate( hzzu);
  deallocate(  kxx);  deallocate(  kxy);  deallocate(  kxz);  deallocate(  kyy);
  deallocate(  kyz);  deallocate(  kzz);  !deallocate(rho_arr);

END PROGRAM coc2cac
