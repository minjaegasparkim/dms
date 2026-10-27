subroutine pre_cals

  !---------------------------------------------------------------------------------------------------|
  ! Pre-calculations                                                                                  |
  !   1. observing wavelengths                                                                        |
  !   2. grain-size grid: n_grain_radius logarithmic bins between grain_radius_min and               |
  !      grain_radius_max; representative radius = geometric bin centre; number of grains in a       |
  !      bin = integral of a^(-q) over the bin (exact), i.e. the size distribution is integrated     |
  !      instead of the former Riemann sum rad^gamma * delrad (whose last bin reached beyond a_max). |
  !   3. single-grain masses                                                                          |
  !---------------------------------------------------------------------------------------------------|

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  integer  :: il, ig, ic, nrow
  real(r2) :: lmin, lmax, step, a1, a2
  real(r2), allocatable :: tab(:,:)

  if (verbose) print *, " "
  if (verbose) print *, " 1) Pre-calculations"
  if (verbose) print *, " "

  !------------------------------------------------------------------------------------------------
  ! 1. observing wavelengths
  !------------------------------------------------------------------------------------------------

  if (switch_wavelength_obs == 2) then
    if (n_wavelength_obs == 1) then
      wavelength_obs(1) = wavelength_obs_min
    else
      lmin = log10(wavelength_obs_min)
      lmax = log10(wavelength_obs_max)
      step = (lmax - lmin) / real(n_wavelength_obs - 1, kind=r2)
      do il = 1, n_wavelength_obs
        wavelength_obs(il) = 10.0_r2**(lmin + real(il - 1, kind=r2) * step)
      end do
    end if
  end if

  allocate(lambda_obs_m(n_wavelength_obs))
  lambda_obs_m(:) = wavelength_obs(:) * micron

  if (verbose) print *, "    observing wavelengths [micron]:"
  do il = 1, n_wavelength_obs
    if (verbose) print *, "      ", il, real_to_string(wavelength_obs(il))
  end do

  !------------------------------------------------------------------------------------------------
  ! 2. grain sizes
  !------------------------------------------------------------------------------------------------

  if (switch_singlegrain_or_Dist == 3) then
    ! size distribution from a file: radius [micron], relative number of grains in the bin
    call read_table(trim(sizedist_file), 2, tab, nrow, min_rows=1)
    if (dms_ierr /= 0) return
    n_grain_radius = nrow
  end if

  allocate(rad(n_grain_radius), rad_edges(0:n_grain_radius), size_weight(n_grain_radius))

  if (switch_singlegrain_or_Dist == 3) then

    rad(:)         = tab(1, :)
    size_weight(:) = tab(2, :)
    do ig = 2, n_grain_radius
      if (rad(ig) <= rad(ig-1)) then
        call dms_fail("the grain radii in "//trim(sizedist_file)//" must increase strictly");  return
      end if
    end do
    if (any(rad <= 0.0_r2) .or. any(size_weight < 0.0_r2) .or. sum(size_weight) <= 0.0_r2) then
      call dms_fail("invalid radii or weights in "//trim(sizedist_file));  return
    end if
    ! bin edges: geometric means of neighbouring radii (first/last bin symmetric in log a)
    do ig = 1, n_grain_radius - 1
      rad_edges(ig) = sqrt(rad(ig) * rad(ig+1))
    end do
    if (n_grain_radius > 1) then
      rad_edges(0)              = rad(1)**2 / rad_edges(1)
      rad_edges(n_grain_radius) = rad(n_grain_radius)**2 / rad_edges(n_grain_radius-1)
    else
      rad_edges(0) = rad(1);  rad_edges(1) = rad(1)
    end if
    size_weight(:)   = size_weight(:) / sum(size_weight)
    grain_radius_min = rad_edges(0)
    grain_radius_max = rad_edges(n_grain_radius)
    q_size           = 0.0_r2

  else if (switch_singlegrain_or_Dist == 1) then

    rad(1)          = grain_radius_min
    rad_edges(0)    = grain_radius_min
    rad_edges(1)    = grain_radius_min
    size_weight(1)  = 1.0_r2

  else

    lmin = log10(grain_radius_min)
    lmax = log10(grain_radius_max)
    do ig = 0, n_grain_radius
      rad_edges(ig) = 10.0_r2**(lmin + (lmax - lmin) * real(ig, kind=r2) / real(n_grain_radius, kind=r2))
    end do

    do ig = 1, n_grain_radius
      a1 = rad_edges(ig - 1)
      a2 = rad_edges(ig)
      rad(ig) = sqrt(a1 * a2)
      if (abs(q_size - 1.0_r2) < 1.0e-10_r2) then
        size_weight(ig) = log(a2 / a1)
      else
        size_weight(ig) = (a2**(1.0_r2 - q_size) - a1**(1.0_r2 - q_size)) / (1.0_r2 - q_size)
      end if
    end do

    size_weight(:) = size_weight(:) / sum(size_weight)

  end if

  if (verbose) print *, " "
  if (verbose) print *, "    grain radii [micron]: ", n_grain_radius, " bin(s) from ", &
           real_to_string(rad(1)), " to ", real_to_string(rad(n_grain_radius))

  !------------------------------------------------------------------------------------------------
  ! 3. grain masses [kg]
  !------------------------------------------------------------------------------------------------

  allocate(grain_mass(n_grain_radius, n_grain_components))

  do ic = 1, n_grain_components
    grain_mass(:, ic) = (4.0_r2 / 3.0_r2) * pi * (rad(:) * micron)**3 &
                        * particle_bulk_density(ic) * 1000.0_r2
  end do

end subroutine pre_cals
