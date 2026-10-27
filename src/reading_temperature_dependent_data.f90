subroutine reading_temperature_dependent_data

  !---------------------------------------------------------------------------------------------------
  ! Reads the optical constants of all dust components and sets up the "optical species".
  !
  ! A component normally has one species (its ri-data file). If temperature-dependent optical
  ! constants are used for component Tdep_comp, that component gets n_Tdep species, one for each
  ! data set (temperature tag + file). A grain then uses, at every distance from the star, the data
  ! set whose temperature is closest to its own equilibrium temperature (self-consistently, also in
  ! the temperature calculation). This replaces the former hard-coded option 1 (crystalline ice at
  ! 15, 25, ..., 170 K, Reinert et al. 2015) and the hard-coded file names of option 2, which made
  ! up about 15000 repetitive lines in the previous version of this file.
  !
  ! The wavelength grid of the optical data ("ri grid") is the union of the wavelengths of all
  ! files, restricted to the range covered by all files; every data set is interpolated onto it
  ! (linear in log lambda). (Up to 2026-10-25 it was the grid of the first file only, so that a
  ! finely sampled band of another component, e.g. the 3 micron ice band in an external mixture
  ! listed after a coarsely sampled refractory material, was smoothed.)
  ! The files may have any number of lines; headers and comment lines are skipped.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  type table_t
    real(r2), allocatable, dimension(:,:) :: t
    integer :: n
  end type table_t

  type(table_t), allocatable, dimension(:) :: tabs
  integer  :: ic, is, k, kk, ia, ilr, n0, i1, i2
  real(r2) :: lam_lo, lam_hi, lam_new
  real(r2), allocatable, dimension(:) :: grid_all, grid_new

  if (verbose) print *, " "
  if (verbose) print *, " 2) Reading optical data"
  if (verbose) print *, " "

  !------------------------------------------------------------------------------------------------
  ! species list
  !------------------------------------------------------------------------------------------------

  allocate(comp_first_species(n_grain_components), comp_n_species(n_grain_components))

  n_species = 0
  do ic = 1, n_grain_components
    comp_first_species(ic) = n_species + 1
    if (switch_Tdep_optics == 1 .and. ic == Tdep_comp) then
      comp_n_species(ic) = n_Tdep
    else
      comp_n_species(ic) = 1
    end if
    n_species = n_species + comp_n_species(ic)
  end do

  allocate(species_comp(n_species), species_T(n_species), species_file(n_species))

  do ic = 1, n_grain_components
    do k = 1, comp_n_species(ic)
      is = comp_first_species(ic) + k - 1
      species_comp(is) = ic
      if (switch_Tdep_optics == 1 .and. ic == Tdep_comp) then
        species_T(is)    = Tdep_T(k)
        species_file(is) = Tdep_file(k)
      else
        species_T(is)    = -1.0_r2
        species_file(is) = dustfile_name(ic)
      end if
    end do
  end do

  !------------------------------------------------------------------------------------------------
  ! read all files
  !------------------------------------------------------------------------------------------------

  allocate(tabs(n_species))
  do is = 1, n_species
    call read_table(trim(ri_dir)//trim(species_file(is)), 3, tabs(is)%t, tabs(is)%n)
    if (dms_ierr /= 0) return
    if (tabs(is)%t(1, tabs(is)%n) < tabs(is)%t(1, 1)) then          ! make ascending
      tabs(is)%t = tabs(is)%t(:, tabs(is)%n:1:-1)
    end if
    if (verbose) print *, "    species ", is, ": component ", species_comp(is), "  ", trim(species_file(is)), &
             "  (", tabs(is)%n, " wavelengths)"
    if (species_T(is) > 0.0_r2 .and. verbose) print *, "                 temperature tag ", species_T(is), " K"
  end do

  lam_lo = -huge(1.0_r2)
  lam_hi =  huge(1.0_r2)
  do is = 1, n_species
    lam_lo = max(lam_lo, tabs(is)%t(1, 1))
    lam_hi = min(lam_hi, tabs(is)%t(1, tabs(is)%n))
  end do

  if (lam_hi <= lam_lo) then
    call dms_fail("the optical-data files have no common wavelength range")
    return
  end if

  ! ri grid = union of the wavelengths of all files within [lam_lo, lam_hi]
  n0 = 0
  do is = 1, n_species
    n0 = n0 + tabs(is)%n
  end do
  allocate(grid_all(n0), grid_new(n0))
  n_wavelength_int = 0
  do is = 1, n_species
    ! wavelengths of this file inside the common range (ascending)
    i1 = 1
    do while (i1 <= tabs(is)%n)
      if (tabs(is)%t(1, i1) >= lam_lo * (1.0_r2 - 1.0e-12_r2)) exit
      i1 = i1 + 1
    end do
    i2 = tabs(is)%n
    do while (i2 >= 1)
      if (tabs(is)%t(1, i2) <= lam_hi * (1.0_r2 + 1.0e-12_r2)) exit
      i2 = i2 - 1
    end do
    if (i2 < i1) cycle
    ! merge with the grid so far, dropping duplicates (relative difference < 1e-9)
    ia = 1;  k = i1;  kk = 0
    do while (ia <= n_wavelength_int .or. k <= i2)
      if (k > i2) then
        lam_new = grid_all(ia);  ia = ia + 1
      else if (ia > n_wavelength_int) then
        lam_new = tabs(is)%t(1, k);  k = k + 1
      else if (grid_all(ia) <= tabs(is)%t(1, k)) then
        lam_new = grid_all(ia);  ia = ia + 1
      else
        lam_new = tabs(is)%t(1, k);  k = k + 1
      end if
      if (kk > 0) then
        if (abs(lam_new - grid_new(kk)) <= 1.0e-9_r2 * lam_new) cycle
      end if
      kk = kk + 1
      grid_new(kk) = lam_new
    end do
    n_wavelength_int = kk
    grid_all(1:kk) = grid_new(1:kk)
  end do

  allocate(wavelength_int(n_wavelength_int), dlambda_int(n_wavelength_int))
  allocate(n_ri(n_wavelength_int, n_species), k_ri(n_wavelength_int, n_species))

  wavelength_int(:) = grid_all(1:n_wavelength_int)
  deallocate(grid_all, grid_new)

  do is = 1, n_species
    do ilr = 1, n_wavelength_int
      n_ri(ilr, is) = interp_lin(log(tabs(is)%t(1,:)), tabs(is)%t(2,:), log(wavelength_int(ilr)))
      k_ri(ilr, is) = interp_lin(log(tabs(is)%t(1,:)), tabs(is)%t(3,:), log(wavelength_int(ilr)))
    end do
  end do

  ! trapezoidal integration weights [m]
  do ilr = 1, n_wavelength_int
    if (n_wavelength_int == 1) then
      dlambda_int(ilr) = 0.0_r2
    else if (ilr == 1) then
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(2) - wavelength_int(1)) * micron
    else if (ilr == n_wavelength_int) then
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(ilr) - wavelength_int(ilr-1)) * micron
    else
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(ilr+1) - wavelength_int(ilr-1)) * micron
    end if
  end do

  if (verbose) print *, " "
  if (verbose) print *, "    common wavelength grid: ", n_wavelength_int, " points, ", &
           real_to_string(wavelength_int(1)), " - ", real_to_string(wavelength_int(n_wavelength_int)), &
           " micron"

  !------------------------------------------------------------------------------------------------
  ! checks
  !------------------------------------------------------------------------------------------------

  if (minval(wavelength_obs) < wavelength_int(1) * (1.0_r2 - 1.0e-9_r2) .or. &
      maxval(wavelength_obs) > wavelength_int(n_wavelength_int) * (1.0_r2 + 1.0e-9_r2)) then
    call dms_fail("the observing wavelengths must lie within the wavelength range of the optical data ("// &
                  real_to_string(wavelength_int(1))//" - "//real_to_string(wavelength_int(n_wavelength_int))// &
                  " micron)")
    return
  end if

  if (stellar_radiance(wavelength_int(1) * micron) > 0.0_r2 .and. &
      wavelength_int(1) > 0.3_r2 * 2897.77_r2 / max(Tstar, 1.0_r2)) then
    if (verbose) print *, "    warning: the optical data start at ", real_to_string(wavelength_int(1)), &
             " micron; part of the stellar UV/optical flux is not absorbed by the grains."
  end if

  deallocate(tabs)

end subroutine reading_temperature_dependent_data
