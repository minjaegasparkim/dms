!____________________________________________________________________________________________________
!
! Input set-up shared by the interactive questions (input_questions) and the input file
! (read_input_file):
!   set_defaults     default values of all input parameters
!   finalize_input   unit conversion, derived quantities and checks
!
! Convention: the input routines store lengths in au and angles in degrees in the belt arrays
! (belt_Rmin, belt_Rmax, belt_Rc, belt_Rstart, belt_Rpeak, belt_mu, belt_sigma, belt_omega);
! finalize_input converts them to m and rad. finalize_input is called exactly once per run.
!____________________________________________________________________________________________________

subroutine set_defaults

  use datatype
  use constants
  use parameters

  implicit none

  dms_ierr = 0;  dms_errmsg = "";  run_complete = .false.

  ri_dir      = "data/optical_constants/"
  results_dir = "results/"
  miex_dir    = "cache/"

  switch_image_or_SED = 1;  switch_reem_or_scat = 3

  switch_BB_or_readingspec = 1
  Tstar = 5800.0_r2;  Rstar_Rsun = 1.0_r2;  Mstar_Msun = 1.0_r2;  dist_star_pc = 10.0_r2
  spec_star = "";  stellar_grid_dir = "";  logg_star = 4.5_r2

  n_belts = 1
  belt_law = 1;  belt_alpha = 1.0_r2;  belt_alpha_out = -3.0_r2
  belt_Rc = 0.0_r2;  belt_Rstart = 0.0_r2;  belt_Rpeak = 0.0_r2;  belt_mu = 0.0_r2;  belt_sigma = 0.0_r2
  belt_Rmin = 0.0_r2;  belt_Rmax = 0.0_r2;  belt_ecc = 0.0_r2;  belt_omega = 0.0_r2
  belt_mass_frac = 1.0_r2;  belt_amp = 1.0_r2
  user_par = 0.0_r2;  user_density_symmetric = .false.

  switch_vertical = 1;  open_angle_deg = 5.0_r2;  h_over_r = 0.05_r2;  gauss_cut = 4.0_r2
  disk_mass_Msun = 1.0e-8_r2;  inc_deg = 0.0_r2;  pa_deg = 0.0_r2

  number_of_pixels = 101;  maxdelT = 2.0_r2;  Nmax = 3;  switch_add_star = 0;  save_cubes = .false.
  switch_psf = 0;  tel_diameter = 0.0_r2;  mask_radius_arcsec = 0.0_r2;  switch_image_unit = 1

  switch_wavelength_obs = 1;  n_wavelength_obs = 0
  wavelength_obs_min = 0.0_r2;  wavelength_obs_max = 0.0_r2

  switch_dust_mix = 0;  switch_dust_mix_method = 1;  vol_frac = 0.0_r2
  emt_matrix_file = "";  emt_inclusion_file = "";  emt_output_file = ""

  n_grain_components = 0
  switch_Tdep_optics = 0;  Tdep_comp = 0;  n_Tdep = 0
  switch_sublimation = 0;  switch_blowout = 1
  switch_blowout_release = 1;  blowout_beta_user = -1.0_r2
  switch_grain_shape = 0;  dhs_fmax = 0.8_r2;  dhs_nf = 10

  switch_singlegrain_or_Dist = 2
  grain_radius_min = 1.0_r2;  grain_radius_max = 1000.0_r2;  q_size = 3.5_r2;  n_grain_radius = 100
  sizedist_file = ""

  n_angle2 = 181
  switch_phase_function = 1;  hg_g1 = 0.0_r2;  hg_g2 = 0.0_r2;  hg_w = 1.0_r2;  pol_max = 0.0_r2

  switch_precal = 0;  switch_miex_save = 0;  file_saved_miex_name = ""
  switch_temperature_distribution = 1;  switch_grainsize_distribution = 1
  switch_output_format = 2;  switch_file_extensions = 0;  file_extension_name = ""

end subroutine set_defaults


subroutine finalize_input

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  integer  :: k, ic
  real(r2) :: fsum

  do_reem = (switch_reem_or_scat == 1 .or. switch_reem_or_scat == 3)
  do_sca  = (switch_reem_or_scat >= 2)

  ! directories end with "/"
  call add_slash(ri_dir);  call add_slash(results_dir);  call add_slash(miex_dir)

  !------------------------------------------------------------------------------------------------
  ! star
  !------------------------------------------------------------------------------------------------
  if (Rstar_Rsun <= 0.0_r2)   call fail_return("the stellar radius must be positive")
  if (dist_star_pc <= 0.0_r2) call fail_return("the distance must be positive")
  if (switch_BB_or_readingspec == 1 .and. Tstar <= 0.0_r2) call fail_return("T_eff must be positive")
  if (switch_BB_or_readingspec == 3 .and. (Tstar <= 0.0_r2 .or. len_trim(stellar_grid_dir) == 0)) &
    call fail_return("a model-atmosphere grid needs T_eff and the grid directory")
  if (dms_ierr /= 0) return
  Rstar     = Rstar_Rsun * R_sun
  Mstar     = Mstar_Msun * M_sun
  dist_star = dist_star_pc * pc

  !------------------------------------------------------------------------------------------------
  ! belts
  !------------------------------------------------------------------------------------------------
  if (n_belts < 1 .or. n_belts > max_belts) then
    call fail_return("the number of belts must be between 1 and "//integer_to_string(max_belts));  return
  end if

  do k = 1, n_belts
    if (belt_law(k) < 1 .or. belt_law(k) > 5) call fail_return("belt "//integer_to_string(k)//": unknown density law")
    if (belt_Rmax(k) <= belt_Rmin(k) .or. belt_Rmin(k) <= 0.0_r2) &
      call fail_return("belt "//integer_to_string(k)//": need 0 < r_in < r_out")
    if (belt_ecc(k) < 0.0_r2 .or. belt_ecc(k) >= 1.0_r2) &
      call fail_return("belt "//integer_to_string(k)//": the eccentricity must be in [0, 1[")
    if (belt_law(k) == 2 .and. belt_Rpeak(k) <= belt_Rstart(k)) &
      call fail_return("belt "//integer_to_string(k)//": R_peak must be larger than R_start")
    if (belt_law(k) == 3 .and. belt_sigma(k) <= 0.0_r2) &
      call fail_return("belt "//integer_to_string(k)//": sigma must be positive")
    if (belt_law(k) == 4 .and. belt_Rc(k) <= 0.0_r2) &
      call fail_return("belt "//integer_to_string(k)//": R_c must be positive")
    if (belt_mass_frac(k) <= 0.0_r2) &
      call fail_return("belt "//integer_to_string(k)//": the mass fraction must be positive")
    if (dms_ierr /= 0) return

    belt_Rmin(k)   = belt_Rmin(k)   * au
    belt_Rmax(k)   = belt_Rmax(k)   * au
    belt_Rc(k)     = belt_Rc(k)     * au
    belt_Rstart(k) = belt_Rstart(k) * au
    belt_Rpeak(k)  = belt_Rpeak(k)  * au
    belt_mu(k)     = belt_mu(k)     * au
    belt_sigma(k)  = belt_sigma(k)  * au
    belt_omega(k)  = belt_omega(k)  * deg2rad
    belt_cos_om(k) = cos(belt_omega(k));  if (abs(belt_cos_om(k)) < 1.0e-14_r2) belt_cos_om(k) = 0.0_r2
    belt_sin_om(k) = sin(belt_omega(k));  if (abs(belt_sin_om(k)) < 1.0e-14_r2) belt_sin_om(k) = 0.0_r2
    belt_rin_eff(k)  = belt_Rmin(k) * (1.0_r2 - belt_ecc(k))
    belt_rout_eff(k) = belt_Rmax(k) * (1.0_r2 + belt_ecc(k))
    ! smallest radial width of the belt (at pericentre), used to sample narrow rings
    belt_width(k) = (belt_Rmax(k) - belt_Rmin(k)) * (1.0_r2 - belt_ecc(k))
    if (belt_law(k) == 3) belt_width(k) = min(belt_width(k), 2.0_r2 * belt_sigma(k) * (1.0_r2 - belt_ecc(k)))
    if (belt_law(k) == 2) belt_width(k) = min(belt_width(k), (belt_Rpeak(k) - belt_Rstart(k)) * (1.0_r2 - belt_ecc(k)))
  end do
  fsum = sum(belt_mass_frac(1:n_belts))
  belt_mass_frac(1:n_belts) = belt_mass_frac(1:n_belts) / fsum
  belt_amp(:) = 1.0_r2

  ! copies of belt 1 (file names, summary)
  switch_density_fn = belt_law(1)
  alpha = belt_alpha(1);  alpha_out = belt_alpha_out(1)
  Rmin = belt_Rmin(1);  Rmax = belt_Rmax(1);  Rmin_au = Rmin / au;  Rmax_au = Rmax / au
  ecc = belt_ecc(1);  omega_peri = belt_omega(1);  omega_peri_deg = omega_peri / deg2rad

  !------------------------------------------------------------------------------------------------
  ! vertical structure, orientation
  !------------------------------------------------------------------------------------------------
  if (switch_vertical == 1) then
    if (open_angle_deg <= 0.0_r2 .or. open_angle_deg >= 90.0_r2) &
      call fail_return("the half opening angle must be in ]0, 90[ deg")
  else if (switch_vertical == 2) then
    if (h_over_r <= 0.0_r2 .or. gauss_cut <= 0.0_r2) call fail_return("h/r and the cut must be positive")
    open_angle_deg = atan(gauss_cut * h_over_r) / deg2rad
  else
    call fail_return("unknown vertical structure")
  end if
  if (dms_ierr /= 0) return

  open_angle = open_angle_deg * deg2rad
  if (.not. (disk_mass_Msun > 0.0_r2)) then
    call dms_fail("the dust mass must be positive");  return
  end if
  disk_mass  = disk_mass_Msun * M_sun
  inc = inc_deg * deg2rad
  pa  = pa_deg  * deg2rad

  sin_open = sin(open_angle)
  tan_open = tan(open_angle)
  cos_inc  = cos(inc);        sin_inc = sin(inc)
  cos_pa   = cos(pa);         sin_pa  = sin(pa)
  cos_om   = belt_cos_om(1);  sin_om  = belt_sin_om(1)
  if (abs(cos_inc) < 1.0e-14_r2) cos_inc = 0.0_r2
  if (abs(sin_inc) < 1.0e-14_r2) sin_inc = 0.0_r2
  if (abs(sin_pa)  < 1.0e-14_r2) sin_pa  = 0.0_r2
  if (abs(cos_pa)  < 1.0e-14_r2) cos_pa  = 0.0_r2

  r_in_eff  = minval(belt_rin_eff(1:n_belts))
  r_out_eff = maxval(belt_rout_eff(1:n_belts))
  if (switch_vertical == 1) then
    r_max_3d = r_out_eff
  else
    r_max_3d = r_out_eff / cos(open_angle)
  end if
  R_field = r_max_3d

  !------------------------------------------------------------------------------------------------
  ! image, wavelengths, dust
  !------------------------------------------------------------------------------------------------
  if (mod(number_of_pixels, 2) == 0) number_of_pixels = number_of_pixels + 1
  if (number_of_pixels < 3) number_of_pixels = 3
  if (Nmax < 0 .or. Nmax > 30) call fail_return("the maximum number of splittings must be in 0 ... 30")
  if (switch_psf == 2 .and. tel_diameter <= 0.0_r2) call fail_return("diffraction PSF: telescope diameter missing")

  if (n_wavelength_obs < 1) call fail_return("no observing wavelengths given")
  if (switch_wavelength_obs == 2 .and. (wavelength_obs_min <= 0.0_r2 .or. wavelength_obs_max < wavelength_obs_min)) &
    call fail_return("invalid wavelength range")
  if (n_grain_components < 1) call fail_return("no dust component given")
  if (dms_ierr /= 0) return

  if (switch_dust_mix == 1 .and. (vol_frac <= 0.0_r2 .or. vol_frac >= 1.0_r2)) &
    call fail_return("EMT: the volume fraction must be in ]0, 1[")
  if (sum(abun(1:n_grain_components)) <= 0.0_r2) call fail_return("the abundances sum up to zero")
  if (any(particle_bulk_density(1:n_grain_components) <= 0.0_r2)) call fail_return("bulk densities must be positive")
  if (switch_Tdep_optics == 1 .and. (Tdep_comp < 1 .or. Tdep_comp > n_grain_components .or. n_Tdep < 1)) &
    call fail_return("temperature-dependent optical data: invalid component or no data sets")
  if (switch_singlegrain_or_Dist == 2 .and. (grain_radius_min <= 0.0_r2 .or. grain_radius_max <= grain_radius_min &
      .or. n_grain_radius < 1)) call fail_return("invalid grain size distribution")
  if (switch_singlegrain_or_Dist == 1 .and. grain_radius_min <= 0.0_r2) call fail_return("invalid grain radius")
  if (switch_singlegrain_or_Dist == 3 .and. len_trim(sizedist_file) == 0) call fail_return("no size-distribution file")
  if (abs(hg_g1) >= 1.0_r2 .or. abs(hg_g2) >= 1.0_r2) call fail_return("|g| must be < 1")
  if (dms_ierr /= 0) return
  abun(1:n_grain_components) = abun(1:n_grain_components) / sum(abun(1:n_grain_components))

  if (switch_sublimation == 1) then
    do ic = 1, n_grain_components
      if (subl_to_comp(ic) < 0 .or. subl_to_comp(ic) > n_grain_components) &
        call fail_return("invalid sublimation residue of component "//integer_to_string(ic))
    end do
    if (dms_ierr /= 0) return
  end if

  if (mod(n_angle2, 2) == 0) n_angle2 = n_angle2 + 1
  if (n_angle2 < 3) n_angle2 = 3
  n_angle = (n_angle2 - 1) / 2 + 1

  if ((switch_precal == 1 .or. switch_miex_save == 1) .and. len_trim(file_saved_miex_name) == 0) &
    file_saved_miex_name = "dms_optics"

  if (switch_image_or_SED == 2) then
    switch_add_star = 0;  maxdelT = 1.0e30_r2;  Nmax = 0;  save_cubes = .false.
  end if

contains

  subroutine fail_return(msg)
    character(len=*), intent(in) :: msg
    call dms_fail(msg)
  end subroutine fail_return

  subroutine add_slash(d)
    character(len=*), intent(inout) :: d
    integer :: n
    n = len_trim(d)
    if (n == 0) then
      d = "./"
    else if (d(n:n) /= "/") then
      d = d(1:n)//"/"
    end if
  end subroutine add_slash

end subroutine finalize_input
