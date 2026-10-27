subroutine max_temp_diff_cals

  !----------------------------------------------------------------------------------------------------
  ! Information only: temperature difference of the smallest grains of component 1 across the
  ! finest (fully split) volume element at the pericentre of the inner edge, i.e. the largest
  ! temperature difference that can remain inside one sub-volume.
  !
  !   maxTempDiff = T(r_in) - T(r_in + sqrt(3) * pixel size / 2^Nmax)
  !----------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2) :: diag

  pixel_size_m = 2.0_r2 * R_field / real(number_of_pixels, kind=r2)
  diag         = sqrt(3.0_r2) * pixel_size_m / 2.0_r2**Nmax

  maxTempDiff = T1_dist(dist_index(r_in_eff)) - T1_dist(dist_index(r_in_eff + diag))

  if (verbose) print *, " "
  if (verbose) print *, "    pixel size: ", real_to_string(pixel_size_m / au), " au;  finest sub-volume: ", &
           real_to_string(pixel_size_m / 2.0_r2**Nmax / au), " au"
  if (verbose) print *, "    max. temperature difference inside the finest sub-volume at the inner edge: ", &
           real_to_string(maxTempDiff), " K"
  if (maxTempDiff > maxdelT) then
    if (verbose) print *, "    note: larger than the requested maxdelT = ", real_to_string(maxdelT), &
             " K; increase Nmax or the number of pixels for a finer sampling."
  end if
  if (verbose) print *, " "

end subroutine max_temp_diff_cals
