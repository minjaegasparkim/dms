subroutine input_questions

  !---------------------------------------------------------------------------------------|
  ! Interactive input of DMS (used when dms is started without an input file).            |
  ! All answers are read with list-directed READ from standard input, so the answers can  |
  ! also be piped from a file:   dms < answers.txt                                        |
  ! The recommended way is an input file (namelist):   dms model.nml                      |
  ! (see examples/ and docs/input_reference.md).                                          |
  !---------------------------------------------------------------------------------------|

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  integer :: ic, k
  logical :: ex

  call set_defaults

  print *, " "
  print *, "============================================================================"
  print *, " "
  print *, "               DMS : Debris disks around Main-sequence Stars                "
  print *, " "
  print *, "   (tip: an input file is easier to reuse:  dms examples/debris_ring.nml)    "
  print *, "============================================================================"
  print *, " "

  !--------------------------------------------------------------------------------
  ! 1. What to compute
  !--------------------------------------------------------------------------------

  print *, "[1] Compute images <1>, SEDs <2>, or both images and SEDs <3> ?"
  read *, switch_image_or_SED
  call check_int(switch_image_or_SED, 1, 3, "image/SED switch");  if (dms_ierr /= 0) return

  print *, "[2] Compute <1> only re-emission, <2> only scattering, <3> both,"
  print *, "    <4> only (linearly polarised) scattering ?"
  read *, switch_reem_or_scat
  call check_int(switch_reem_or_scat, 1, 4, "re-emission/scattering switch");  if (dms_ierr /= 0) return

  !--------------------------------------------------------------------------------
  ! 2. Star
  !--------------------------------------------------------------------------------

  print *, "[3] Stellar spectrum: <1> black body, <2> spectrum from a file, <3> model-atmosphere grid ?"
  read *, switch_BB_or_readingspec
  call check_int(switch_BB_or_readingspec, 1, 3, "stellar spectrum switch");  if (dms_ierr /= 0) return

  select case (switch_BB_or_readingspec)
  case (1)
    print *, "[3a] Effective temperature of the star [K] ?"
    read *, Tstar
  case (2)
    print *, "[3a] Name of the spectrum file (relative to the working directory). Columns:"
    print *, "     wavelength [micron], flux density [Jy] as observed at the distance of the star."
    read *, spec_star
    inquire(file=trim(spec_star), exist=ex)
    if (.not. ex) then
      call dms_fail("stellar spectrum file not found: "//trim(spec_star));  return
    end if
  case (3)
    print *, "[3a] Directory of the model-atmosphere grid (with grid_index.txt) ?"
    read *, stellar_grid_dir
    print *, "[3b] Effective temperature [K] and log g [cgs] ?"
    read *, Tstar, logg_star
  end select

  print *, "[4] Stellar radius [R_sun] ?"
  read *, Rstar_Rsun

  print *, "[5] Stellar mass [M_sun] ? (used for the blow-out size; any value if not needed)"
  read *, Mstar_Msun

  print *, "[6] Distance of the star [pc] ?"
  read *, dist_star_pc

  !--------------------------------------------------------------------------------
  ! 3. Disk: belts
  !--------------------------------------------------------------------------------

  print *, "[7] Number of belts (rings) with different radial distributions (1 ... ", max_belts, ") ?"
  read *, n_belts
  call check_int(n_belts, 1, max_belts, "number of belts");  if (dms_ierr /= 0) return

  do k = 1, n_belts

    if (n_belts > 1) print *, "    --- belt ", k, " ---"

    print *, "[8] Radial density distribution (a = semi-major axis):"
    print *, "    <1> n ~ a^(-alpha)"
    print *, "    <2> n ~ ((a-R_start)/(R_peak-R_start))^alpha exp[alpha (1-(a-R_start)/(R_peak-R_start))]"
    print *, "    <3> n ~ a^alpha exp[-0.5 ((a - mu)/sigma)^2]"
    print *, "    <4> n ~ [(a/R_c)^(-2 alpha_in) + (a/R_c)^(-2 alpha_out)]^(-1/2)  (alpha_in>0, alpha_out<0)"
    print *, "    <5> user-defined (user_density.f90, parameters below)"
    read *, belt_law(k)
    call check_int(belt_law(k), 1, 5, "density distribution");  if (dms_ierr /= 0) return

    select case (belt_law(k))
    case (1)
      print *, "[8a] Exponent alpha ?"
      read *, belt_alpha(k)
    case (2)
      print *, "[8a] Exponent alpha, R_start [au], R_peak [au] ?"
      read *, belt_alpha(k), belt_Rstart(k), belt_Rpeak(k)
    case (3)
      print *, "[8a] Exponent alpha, mu [au], sigma [au] ?"
      read *, belt_alpha(k), belt_mu(k), belt_sigma(k)
    case (4)
      print *, "[8a] alpha_in (> 0), alpha_out (< 0), R_c [au] ?"
      read *, belt_alpha(k), belt_alpha_out(k), belt_Rc(k)
    case (5)
      print *, "[8a] Up to 10 parameters user_par(1:10) for user_density_fn (end the list with /) ?"
      read *, user_par
    end select

    print *, "[9] Inner and outer radius [au] (for eccentric belts: smallest and largest semi-major axis) ?"
    read *, belt_Rmin(k), belt_Rmax(k)

    print *, "[10] Eccentricity of the dust orbits (0 = circular) ?"
    read *, belt_ecc(k)
    belt_omega(k) = 0.0_r2
    if (belt_ecc(k) > 0.0_r2) then
      print *, "[10a] Argument of pericentre omega [deg], measured in the disk plane from the"
      print *, "      line of nodes (+Y, the projected major axis) towards +X, the far side of"
      print *, "      the disk (90 = pericentre on the far side of the minor axis) ?"
      read *, belt_omega(k)
    end if

    belt_mass_frac(k) = 1.0_r2
    if (n_belts > 1) then
      print *, "[11] Relative dust mass of this belt (any unit; normalised over all belts) ?"
      read *, belt_mass_frac(k)
    end if

  end do

  !--------------------------------------------------------------------------------
  ! 4. Disk: vertical structure, mass, orientation
  !--------------------------------------------------------------------------------

  print *, "[12] Vertical structure: <1> wedge (constant density within the opening angle),"
  print *, "     <2> Gaussian with scale height h = (h/r) r ?"
  read *, switch_vertical
  call check_int(switch_vertical, 1, 2, "vertical structure");  if (dms_ierr /= 0) return
  if (switch_vertical == 1) then
    print *, "[12a] Half opening angle of the disk [deg] ?"
    read *, open_angle_deg
  else
    print *, "[12a] Aspect ratio h/r ?"
    read *, h_over_r
  end if

  print *, "[13] Total dust mass [M_sun] (mass of all grains between the minimum and"
  print *, "     maximum grain size; after sublimation/blow-out if these are switched on) ?"
  read *, disk_mass_Msun

  print *, "[14] Inclination of the disk [deg] (0 = face-on, 90 = edge-on) ?"
  read *, inc_deg

  print *, "[15] Position angle of the disk major axis [deg], counter-clockwise from the"
  print *, "     +y axis of the image (0 = major axis along y, as in previous DMS versions) ?"
  read *, pa_deg

  !--------------------------------------------------------------------------------
  ! 5. Images
  !--------------------------------------------------------------------------------

  if (switch_image_or_SED == 1 .or. switch_image_or_SED == 3) then

    print *, "[16] Number of pixels in one direction (odd number) ?"
    read *, number_of_pixels

    print *, "[17] Maximum temperature difference in one (sub-)volume [K] ?"
    read *, maxdelT

    print *, "[18] Maximum number of splittings of one volume element (3-4 is usually enough) ?"
    read *, Nmax
    call check_int(Nmax, 0, 30, "maximum number of splittings");  if (dms_ierr /= 0) return

    print *, "[19] Add the stellar flux to the central pixel of the total image? [1] yes, [0] no"
    read *, switch_add_star

    print *, "[19a] Save 3-d cubes (flux of every volume element and its scattering angle)? [1] yes, [0] no"
    read *, k
    save_cubes = (k == 1)

    print *, "[20] Simulated observation: PSF <0> none, <1> Gaussian (FWHM), <2> diffraction (telescope) ?"
    read *, switch_psf
    call check_int(switch_psf, 0, 2, "PSF switch");  if (dms_ierr /= 0) return
    if (switch_psf == 1) then
      print *, "[20a] FWHM of the PSF [arcsec] (same at all wavelengths) ?"
      read *, var
    else if (switch_psf == 2) then
      print *, "[20a] Telescope diameter [m] (FWHM = 1.029 lambda/D) ?"
      read *, tel_diameter
    end if
    print *, "[20b] Radius of an inner mask (coronagraph) [arcsec] (0 = none) ?"
    read *, mask_radius_arcsec

    print *, "[21] Unit of the images: <1> Jy/pixel, <2> mJy/arcsec^2, <3> Jy/arcsec^2 ?"
    read *, switch_image_unit
    call check_int(switch_image_unit, 1, 3, "image unit");  if (dms_ierr /= 0) return

  end if

  !--------------------------------------------------------------------------------
  ! 6. Observing wavelengths
  !--------------------------------------------------------------------------------

  print *, "[22] Observing wavelengths: list <1> or logarithmic distribution <2> ?"
  read *, switch_wavelength_obs
  call check_int(switch_wavelength_obs, 1, 2, "wavelength switch");  if (dms_ierr /= 0) return

  print *, "[23] Number of observing wavelengths ?"
  read *, n_wavelength_obs
  call check_int(n_wavelength_obs, 1, 100000, "number of wavelengths");  if (dms_ierr /= 0) return
  allocate(wavelength_obs(n_wavelength_obs), psf_fwhm_arcsec(n_wavelength_obs))
  wavelength_obs = 0.0_r2
  psf_fwhm_arcsec = 0.0_r2
  if (switch_psf == 1) psf_fwhm_arcsec = var

  if (switch_wavelength_obs == 1) then
    do k = 1, n_wavelength_obs
      print *, "[23a] Observing wavelength number ", k, " [micron] ?"
      read *, wavelength_obs(k)
    end do
  else
    print *, "[23a] Smallest and largest observing wavelength [micron] ?"
    read *, wavelength_obs_min, wavelength_obs_max
  end if

  !--------------------------------------------------------------------------------
  ! 7. Effective medium theory (optional): creates a new optical-data file
  !--------------------------------------------------------------------------------

  print *, "[24] Create a mixed material with an effective medium approximation? [1] yes, [0] no"
  print *, "     (the result is written to the optical-data directory and can be used below)"
  read *, switch_dust_mix

  if (switch_dust_mix == 1) then
    print *, "[24a] File of the matrix material ?"
    read *, emt_matrix_file
    print *, "[24b] File of the inclusion material ?"
    read *, emt_inclusion_file
    print *, "[24c] Volume fraction of the inclusion (0 < f < 1) ?"
    read *, vol_frac
    print *, "[24d] Mixing rule: <1> Maxwell-Garnett, <2> Bruggeman, <3> Lichtenecker, <4> Looyenga ?"
    read *, switch_dust_mix_method
    call check_int(switch_dust_mix_method, 1, 4, "mixing rule");  if (dms_ierr /= 0) return
    print *, "[24e] Name of the output file ?"
    read *, emt_output_file
  end if

  !--------------------------------------------------------------------------------
  ! 8. Dust components
  !--------------------------------------------------------------------------------

  print *, "[25] Number of dust components (materials) ?"
  read *, n_grain_components
  call check_int(n_grain_components, 1, 100, "number of components");  if (dms_ierr /= 0) return

  allocate(dustfile_name(n_grain_components), particle_bulk_density(n_grain_components), &
           abun(n_grain_components))

  print *, "[26] Optical-data file of each component (in the optical-data directory;"
  print *, "     columns wavelength [micron], n, k)"
  do ic = 1, n_grain_components
    print *, "     component ", ic, " :"
    read *, dustfile_name(ic)
  end do

  print *, "[27] Bulk density of each component [g/cm^3]"
  do ic = 1, n_grain_components
    print *, "     component ", ic, " :"
    read *, particle_bulk_density(ic)
  end do

  abun(:) = 1.0_r2
  if (n_grain_components > 1) then
    print *, "[28] Relative abundance of each component [%] (same size distribution for all;"
    print *, "     a sublimation residue that only forms by sublimation gets 0)"
    do ic = 1, n_grain_components
      print *, "     component ", ic, " :"
      read *, abun(ic)
    end do
  end if

  !--------------------------------------------------------------------------------
  ! 9. Temperature-dependent optical constants
  !--------------------------------------------------------------------------------

  print *, "[29] Use temperature-dependent optical constants for one component? [1] yes, [0] no"
  read *, switch_Tdep_optics

  if (switch_Tdep_optics == 1) then
    print *, "[29a] Which component ?"
    read *, Tdep_comp
    call check_int(Tdep_comp, 1, n_grain_components, "component index");  if (dms_ierr /= 0) return
    print *, "[29b] Number of data sets (temperatures) ?"
    read *, n_Tdep
    call check_int(n_Tdep, 1, 1000, "number of data sets");  if (dms_ierr /= 0) return
    allocate(Tdep_T(n_Tdep), Tdep_file(n_Tdep))
    do k = 1, n_Tdep
      print *, "[29c] Data set ", k, " : temperature [K] and file name ?"
      read *, Tdep_T(k), Tdep_file(k)
    end do
  else
    n_Tdep = 0
  end if

  !--------------------------------------------------------------------------------
  ! 10. Sublimation and radiation-pressure blow-out
  !--------------------------------------------------------------------------------

  print *, "[30] Consider sublimation of the grains? [1] yes, [0] no"
  read *, switch_sublimation

  allocate(T_subl(n_grain_components), subl_to_comp(n_grain_components))
  T_subl(:)       = -1.0_r2
  subl_to_comp(:) = 0

  if (switch_sublimation == 1) then
    print *, "[30a] Sublimation temperature of each component [K] (<= 0: never sublimates)"
    do ic = 1, n_grain_components
      print *, "      component ", ic, " :"
      read *, T_subl(ic)
    end do
    print *, "[30b] What remains after sublimation of each component?"
    print *, "      0: the grain vanishes,  n: the grain becomes component n"
    print *, "      (e.g. ice/silicate mixture -> porous silicate, Kim et al. 2019)"
    do ic = 1, n_grain_components
      print *, "      component ", ic, " :"
      read *, subl_to_comp(ic)
    end do
  end if

  print *, "[31] Remove grains blown out by radiation pressure (beta >= 0.5)? [1] yes, [0] no"
  read *, switch_blowout

  !--------------------------------------------------------------------------------
  ! 11. Grain sizes
  !--------------------------------------------------------------------------------

  print *, "[32] Grain sizes: <1> single size, <2> power-law distribution, <3> from a file ?"
  read *, switch_singlegrain_or_Dist
  call check_int(switch_singlegrain_or_Dist, 1, 3, "grain size switch");  if (dms_ierr /= 0) return

  select case (switch_singlegrain_or_Dist)
  case (1)
    print *, "[32a] Grain radius [micron] ?"
    read *, grain_radius_min
    grain_radius_max = grain_radius_min
    n_grain_radius   = 1
    q_size           = 0.0_r2
  case (2)
    print *, "[32a] Minimum and maximum grain radius [micron] ?"
    read *, grain_radius_min, grain_radius_max
    print *, "[32b] Exponent q of the size distribution n(a) ~ a^(-q) (e.g. 3.5) ?"
    read *, q_size
    print *, "[32c] Number of (logarithmic) grain size bins ?"
    read *, n_grain_radius
  case (3)
    print *, "[32a] File with grain radius [micron] and relative number of grains (per bin) ?"
    read *, sizedist_file
  end select

  !--------------------------------------------------------------------------------
  ! 12. Scattering
  !--------------------------------------------------------------------------------

  print *, "[33] Number of scattering angles between 0 and 180 deg (odd, e.g. 181) ?"
  read *, n_angle2

  print *, "[34] Scattering phase function: <1> Mie theory, <2> Henyey-Greenstein,"
  print *, "     <3> two-component Henyey-Greenstein (e.g. Engler et al. 2017) ?"
  read *, switch_phase_function
  call check_int(switch_phase_function, 1, 3, "phase function switch");  if (dms_ierr /= 0) return

  if (switch_phase_function == 2) then
    print *, "[34a] Asymmetry parameter g ?"
    read *, hg_g1
    print *, "[34b] Maximum degree of linear polarisation p_max (Rayleigh-like, 0..1) ?"
    read *, pol_max
  else if (switch_phase_function == 3) then
    print *, "[34a] g1, g2 and weight w of the first component:  w HG(g1) + (1-w) HG(g2) ?"
    read *, hg_g1, hg_g2, hg_w
    print *, "[34b] Maximum degree of linear polarisation p_max (Rayleigh-like, 0..1) ?"
    read *, pol_max
  end if

  !--------------------------------------------------------------------------------
  ! 13. Saving / loading
  !--------------------------------------------------------------------------------

  print *, "[35] Load previously saved optical properties (Mie results) <1> or compute new ones <0> ?"
  read *, switch_precal

  print *, "[36] Save the optical properties for later runs? [1] yes, [0] no"
  read *, switch_miex_save

  if (switch_precal == 1 .or. switch_miex_save == 1) then
    print *, "[36a] Name of the optical-properties file (in the cache directory) ?"
    read *, file_saved_miex_name
  end if

  print *, "[37] Save the temperature distribution? [1] yes, [0] no"
  read *, switch_temperature_distribution

  print *, "[38] Save the grain size distribution (numbers, cross sections, beta)? [1] yes, [0] no"
  read *, switch_grainsize_distribution

  print *, "[39] Output format of images: <1> ASCII, <2> FITS, <3> both ?"
  read *, switch_output_format
  call check_int(switch_output_format, 1, 3, "output format");  if (dms_ierr /= 0) return

  print *, "[40] File extension: [0] none, [1] name + model parameters, [2] name only ?"
  read *, switch_file_extensions
  call check_int(switch_file_extensions, 0, 2, "file extension switch");  if (dms_ierr /= 0) return

  file_extension_name = ""
  if (switch_file_extensions > 0) then
    print *, "[40a] Name (e.g. of your project) ?"
    read *, file_extension_name
  end if

  print *, " "
  print *, "  Input complete - calculation is starting."
  print *, " "

contains

  subroutine check_int(val, lo, hi, what)
    integer, intent(in)          :: val, lo, hi
    character(len=*), intent(in) :: what
    if (val < lo .or. val > hi) then
      call dms_fail("invalid value for the "//what//": "//integer_to_string(val)// &
                    " (allowed "//integer_to_string(lo)//" ... "//integer_to_string(hi)//")")
    end if
  end subroutine check_int

end subroutine input_questions
