subroutine scattering_angles_cals

  !---------------------------------------------------------------------------------------------------
  ! Scattering-angle grid: n_angle2 equidistant angles in [0, 180] deg
  ! (the same grid on which miex returns S1, S2)
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters

  implicit none

  integer :: ia

  allocate(x_angle(n_angle2))

  d_angle = 180.0_r2 / real(n_angle2 - 1, kind=r2)

  do ia = 1, n_angle2
    x_angle(ia) = real(ia - 1, kind=r2) * d_angle
  end do

  if (verbose) print *, " "
  if (verbose) print *, "    scattering angles: ", n_angle2, " values, step ", d_angle, " deg"

end subroutine scattering_angles_cals
