module grid_creation
  use, intrinsic :: iso_fortran_env, only: long => real64
  implicit none
contains

  subroutine create_cartesian_grid(x_cart, y_cart, z_cart, res)
    real(long), pointer :: x_cart(:), y_cart(:), z_cart(:)
    integer, intent(in) :: res
    real(long) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(long) :: dx, dy, dz
    integer :: i

    x_min = -20.0_long
    x_max = 20.0_long
    y_min = -20.0_long
    y_max = 20.0_long
    z_min = -20.0_long
    z_max = 20.0_long

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
  end subroutine create_cartesian_grid

end module grid_creation
