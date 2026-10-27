!____________________________________________________________________________________________________
!
! Input file (Fortran namelist), e.g.   ./dms examples/debris_ring.nml
!____________________________________________________________________________________________________
!
! Groups (all optional, any order; unspecified values keep their defaults):
!   &run          what to compute, directories, output options
!   &star         stellar parameters and spectrum
!   &disk         belts, vertical structure, orientation, dust mass
!   &image        image size, sampling, simulated observations
!   &wavelengths  observing wavelengths
!   &dust         components, optical data, size distribution, EMT, sublimation, blow-out
!   &scattering   angles and phase function
! The complete list with units and defaults is given in docs/input_reference.md.
! Lengths in au, angles in degrees, masses in M_sun, radii of grains in micron.
!____________________________________________________________________________________________________

subroutine read_input_file(fname)

  use datatype
  use constants
  use parameters, pm_alpha => alpha, pm_alpha_out => alpha_out, pm_mu => mu, pm_sigma => sigma, &
                  pm_R_c => R_c, pm_R_start => R_start, pm_R_peak => R_peak, &
                  pm_n_belts => n_belts, pm_h_over_r => h_over_r, pm_gauss_cut => gauss_cut, &
                  pm_user_par => user_par, pm_user_density_symmetric => user_density_symmetric, &
                  pm_dhs_fmax => dhs_fmax, pm_dhs_nf => dhs_nf
  use functions, only: dms_fail, integer_to_string

  implicit none

  character(len=*), intent(in) :: fname

  integer, parameter :: mb = max_belts, mc = 50, mw = 5000, mt = 200

  ! &run
  character(len=32)  :: mode, emission, output_format
  character(len=512) :: data_dir, output_dir, cache_dir, output_name, optics_cache_name
  logical :: name_with_parameters, save_temperature, save_grain_sizes, load_optics, save_optics, &
             screen_output, write_files
  ! &star
  character(len=32)  :: spectrum
  character(len=512) :: spectrum_file, grid_dir
  real(r2) :: teff, radius, mass, distance, logg
  ! &disk
  integer  :: n_belts
  character(len=32), dimension(mb) :: density_law
  real(r2), dimension(mb) :: alpha, alpha_out, r_c, r_start, r_peak, mu, sigma, r_in, r_out, &
                             eccentricity, omega, belt_mass_fraction
  real(r2), dimension(10) :: user_par
  logical  :: user_density_symmetric
  character(len=32) :: vertical
  real(r2) :: opening_angle, h_over_r, gauss_cut, dust_mass, inclination, position_angle
  ! &image
  integer  :: n_pixels, n_split_max
  real(r2) :: max_delta_t, telescope_diameter, mask_radius
  logical  :: add_star, save_3d_cubes
  character(len=32) :: psf, image_unit
  real(r2), dimension(mw) :: psf_fwhm
  ! &wavelengths
  real(r2), dimension(mw) :: lambda
  real(r2) :: lambda_min, lambda_max
  integer  :: n_lambda
  ! &dust
  integer  :: n_components, tdep_component
  character(len=512), dimension(mc) :: optical_data
  real(r2), dimension(mc) :: bulk_density, abundance, t_sublimation
  integer,  dimension(mc) :: sublimation_residue
  character(len=32)  :: size_distribution, emt_rule, grain_shape, blowout_release
  character(len=512) :: size_file, emt_matrix, emt_inclusion, emt_output
  real(r2) :: a_min, a_max, q, emt_volume_fraction, dhs_fmax, blowout_beta
  integer  :: n_sizes, dhs_nf
  logical  :: emt, sublimation, blowout
  real(r2), dimension(mt) :: tdep_temperatures
  character(len=512), dimension(mt) :: tdep_files
  ! &scattering
  integer  :: n_angles
  character(len=32) :: phase_function
  real(r2) :: g1, g2, w, p_max

  namelist /run/ mode, emission, data_dir, output_dir, cache_dir, output_name, output_format, &
                 name_with_parameters, save_temperature, save_grain_sizes, load_optics, save_optics, &
                 optics_cache_name, screen_output, write_files
  namelist /star/ spectrum, teff, radius, mass, distance, logg, spectrum_file, grid_dir
  namelist /disk/ n_belts, density_law, alpha, alpha_out, r_c, r_start, r_peak, mu, sigma, r_in, r_out, &
                  eccentricity, omega, belt_mass_fraction, user_par, user_density_symmetric, vertical, &
                  opening_angle, h_over_r, gauss_cut, dust_mass, inclination, position_angle
  namelist /image/ n_pixels, max_delta_t, n_split_max, add_star, psf, psf_fwhm, telescope_diameter, &
                   mask_radius, image_unit, save_3d_cubes
  namelist /wavelengths/ lambda, lambda_min, lambda_max, n_lambda
  namelist /dust/ n_components, optical_data, bulk_density, abundance, size_distribution, a_min, a_max, q, &
                  n_sizes, size_file, emt, emt_matrix, emt_inclusion, emt_volume_fraction, emt_rule, &
                  emt_output, tdep_component, tdep_temperatures, tdep_files, sublimation, t_sublimation, &
                  sublimation_residue, blowout, blowout_beta, blowout_release, grain_shape, dhs_fmax, dhs_nf
  namelist /scattering/ n_angles, phase_function, g1, g2, w, p_max

  integer :: u, ios, k, nw, nf
  character(len=1024) :: line
  logical :: ex
  character(len=256) :: iomsg_txt

  call set_defaults

  inquire(file=fname, exist=ex)
  if (.not. ex) then
    call dms_fail("input file not found: "//trim(fname));  return
  end if

  !------------------------------------------------------------------------------------------------
  ! defaults of the namelist variables (from set_defaults)
  !------------------------------------------------------------------------------------------------
  mode = "image";  emission = "both";  output_format = "fits"
  data_dir = ri_dir;  output_dir = results_dir;  cache_dir = miex_dir
  output_name = "";  optics_cache_name = "";  name_with_parameters = .false.
  save_temperature = .true.;  save_grain_sizes = .true.;  load_optics = .false.;  save_optics = .false.
  screen_output = .true.;  write_files = .true.

  spectrum = "blackbody";  teff = Tstar;  radius = Rstar_Rsun;  mass = Mstar_Msun;  distance = dist_star_pc
  logg = logg_star;  spectrum_file = "";  grid_dir = ""

  n_belts = 1;  density_law = "powerlaw";  alpha = 1.0_r2;  alpha_out = -3.0_r2
  r_c = 0.0_r2;  r_start = 0.0_r2;  r_peak = 0.0_r2;  mu = 0.0_r2;  sigma = 0.0_r2
  r_in = 0.0_r2;  r_out = 0.0_r2;  eccentricity = 0.0_r2;  omega = 0.0_r2;  belt_mass_fraction = 1.0_r2
  user_par = 0.0_r2;  user_density_symmetric = .false.
  vertical = "wedge";  opening_angle = open_angle_deg;  h_over_r = 0.05_r2;  gauss_cut = 4.0_r2
  dust_mass = disk_mass_Msun;  inclination = 0.0_r2;  position_angle = 0.0_r2

  n_pixels = number_of_pixels;  max_delta_t = maxdelT;  n_split_max = Nmax;  add_star = .false.
  save_3d_cubes = .false.
  psf = "none";  psf_fwhm = -1.0_r2;  telescope_diameter = 0.0_r2;  mask_radius = 0.0_r2
  image_unit = "Jy/pixel"

  lambda = -1.0_r2;  lambda_min = 0.0_r2;  lambda_max = 0.0_r2;  n_lambda = 0

  n_components = 0;  optical_data = "";  bulk_density = 0.0_r2;  abundance = 1.0_r2
  size_distribution = "powerlaw";  a_min = grain_radius_min;  a_max = grain_radius_max;  q = q_size
  n_sizes = n_grain_radius;  size_file = ""
  emt = .false.;  emt_matrix = "";  emt_inclusion = "";  emt_volume_fraction = 0.0_r2
  emt_rule = "maxwell-garnett";  emt_output = ""
  tdep_component = 0;  tdep_temperatures = -1.0_r2;  tdep_files = ""
  sublimation = .false.;  t_sublimation = -1.0_r2;  sublimation_residue = 0;  blowout = .true.
  blowout_beta = -1.0_r2;  blowout_release = 'circular'
  grain_shape = merge('dhs   ', 'sphere', switch_grain_shape == 1);  dhs_fmax = pm_dhs_fmax;  dhs_nf = pm_dhs_nf

  n_angles = n_angle2;  phase_function = "mie";  g1 = 0.0_r2;  g2 = 0.0_r2;  w = 1.0_r2;  p_max = 0.0_r2

  !------------------------------------------------------------------------------------------------
  ! read the groups
  !------------------------------------------------------------------------------------------------
  ! groups with unknown names (e.g. a typo such as &disc) would be ignored silently
  open(newunit=u, file=fname, action="read", status="old", form="formatted")
  do
    read(u, '(a)', iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (line(1:1) /= "&") cycle
    k = scan(line(2:), " /,") 
    if (k == 0) k = len_trim(line)
    select case (lower(line(2:k)))
    case ("run", "star", "disk", "image", "wavelengths", "dust", "scattering")
    case default
      call dms_fail("unknown namelist group "//trim(line(1:k))//" in "//trim(fname)// &
                    " (allowed: &run &star &disk &image &wavelengths &dust &scattering)")
      close(u);  return
    end select
  end do
  close(u)

  open(newunit=u, file=fname, action="read", status="old", form="formatted")
  rewind(u);  read(u, nml=run,         iostat=ios, iomsg=iomsg_txt);  call check("run")
  rewind(u);  read(u, nml=star,        iostat=ios, iomsg=iomsg_txt);  call check("star")
  rewind(u);  read(u, nml=disk,        iostat=ios, iomsg=iomsg_txt);  call check("disk")
  rewind(u);  read(u, nml=image,       iostat=ios, iomsg=iomsg_txt);  call check("image")
  rewind(u);  read(u, nml=wavelengths, iostat=ios, iomsg=iomsg_txt);  call check("wavelengths")
  rewind(u);  read(u, nml=dust,        iostat=ios, iomsg=iomsg_txt);  call check("dust")
  rewind(u);  read(u, nml=scattering,  iostat=ios, iomsg=iomsg_txt);  call check("scattering")
  close(u)
  if (dms_ierr /= 0) return

  !------------------------------------------------------------------------------------------------
  ! &run
  !------------------------------------------------------------------------------------------------
  select case (lower(mode))
  case ("image");       switch_image_or_SED = 1
  case ("sed");         switch_image_or_SED = 2
  case ("both");        switch_image_or_SED = 3
  case default;         call dms_fail("&run mode must be 'image', 'sed' or 'both'");  return
  end select
  select case (lower(emission))
  case ("thermal", "reemission");             switch_reem_or_scat = 1
  case ("scattering");                        switch_reem_or_scat = 2
  case ("both");                              switch_reem_or_scat = 3
  case ("polarisation", "polarization");      switch_reem_or_scat = 4
  case default
    call dms_fail("&run emission must be 'thermal', 'scattering', 'both' or 'polarisation'");  return
  end select
  select case (lower(output_format))
  case ("ascii");  switch_output_format = 1
  case ("fits");   switch_output_format = 2
  case ("both");   switch_output_format = 3
  case default;    call dms_fail("&run output_format must be 'ascii', 'fits' or 'both'");  return
  end select
  ri_dir = data_dir;  results_dir = output_dir;  miex_dir = cache_dir
  file_extension_name = output_name
  switch_file_extensions = 0
  if (len_trim(output_name) > 0) then
    switch_file_extensions = 2
    if (name_with_parameters) switch_file_extensions = 1
  end if
  switch_temperature_distribution = merge(1, 0, save_temperature)
  switch_grainsize_distribution   = merge(1, 0, save_grain_sizes)
  switch_precal    = merge(1, 0, load_optics)
  switch_miex_save = merge(1, 0, save_optics)
  file_saved_miex_name = optics_cache_name
  verbose = screen_output .and. .not. force_quiet
  write_output = write_files

  !------------------------------------------------------------------------------------------------
  ! &star
  !------------------------------------------------------------------------------------------------
  select case (lower(spectrum))
  case ("blackbody", "bb");  switch_BB_or_readingspec = 1
  case ("file");             switch_BB_or_readingspec = 2
  case ("grid");             switch_BB_or_readingspec = 3
  case default;              call dms_fail("&star spectrum must be 'blackbody', 'file' or 'grid'");  return
  end select
  Tstar = teff;  Rstar_Rsun = radius;  Mstar_Msun = mass;  dist_star_pc = distance;  logg_star = logg
  spec_star = spectrum_file;  stellar_grid_dir = grid_dir
  if (switch_BB_or_readingspec == 2) then
    inquire(file=trim(spec_star), exist=ex)
    if (.not. ex) then
      call dms_fail("stellar spectrum file not found: "//trim(spec_star));  return
    end if
  end if

  !------------------------------------------------------------------------------------------------
  ! &disk
  !------------------------------------------------------------------------------------------------
  pm_n_belts = n_belts
  if (n_belts < 1 .or. n_belts > mb) then
    call dms_fail("&disk n_belts must be between 1 and "//integer_to_string(mb));  return
  end if
  do k = 1, n_belts
    select case (lower(density_law(k)))
    case ("powerlaw", "power-law", "1");      belt_law(k) = 1
    case ("ertel", "2");                      belt_law(k) = 2
    case ("gaussian", "3");                   belt_law(k) = 3
    case ("broken", "broken-powerlaw", "4");  belt_law(k) = 4
    case ("user", "5");                       belt_law(k) = 5
    case default
      call dms_fail("&disk density_law("//integer_to_string(k)//") unknown: "//trim(density_law(k)));  return
    end select
    belt_alpha(k) = alpha(k);  belt_alpha_out(k) = alpha_out(k)
    belt_Rc(k) = r_c(k);  belt_Rstart(k) = r_start(k);  belt_Rpeak(k) = r_peak(k)
    belt_mu(k) = mu(k);   belt_sigma(k) = sigma(k)
    belt_Rmin(k) = r_in(k);  belt_Rmax(k) = r_out(k)
    belt_ecc(k) = eccentricity(k);  belt_omega(k) = omega(k)
    belt_mass_frac(k) = belt_mass_fraction(k)
  end do
  pm_user_par = user_par;  pm_user_density_symmetric = user_density_symmetric
  select case (lower(vertical))
  case ("wedge");     switch_vertical = 1
  case ("gaussian");  switch_vertical = 2
  case default;       call dms_fail("&disk vertical must be 'wedge' or 'gaussian'");  return
  end select
  open_angle_deg = opening_angle;  pm_h_over_r = h_over_r;  pm_gauss_cut = gauss_cut
  disk_mass_Msun = dust_mass;  inc_deg = inclination;  pa_deg = position_angle

  !------------------------------------------------------------------------------------------------
  ! &wavelengths (needed before &image for the PSF list)
  !------------------------------------------------------------------------------------------------
  nw = count(lambda > 0.0_r2)
  if (nw > 0) then
    switch_wavelength_obs = 1
    n_wavelength_obs = nw
    allocate(wavelength_obs(nw))
    wavelength_obs(:) = pack(lambda, lambda > 0.0_r2)
  else if (n_lambda > 0) then
    switch_wavelength_obs = 2
    n_wavelength_obs = n_lambda
    wavelength_obs_min = lambda_min;  wavelength_obs_max = lambda_max
    allocate(wavelength_obs(n_lambda))
    wavelength_obs = 0.0_r2
  else
    call dms_fail("&wavelengths: give either lambda = ... (list) or lambda_min, lambda_max, n_lambda");  return
  end if

  !------------------------------------------------------------------------------------------------
  ! &image
  !------------------------------------------------------------------------------------------------
  number_of_pixels = n_pixels;  maxdelT = max_delta_t;  Nmax = n_split_max
  switch_add_star = merge(1, 0, add_star)
  save_cubes = save_3d_cubes
  select case (lower(psf))
  case ("none");         switch_psf = 0
  case ("gaussian");     switch_psf = 1
  case ("diffraction");  switch_psf = 2
  case default;          call dms_fail("&image psf must be 'none', 'gaussian' or 'diffraction'");  return
  end select
  allocate(psf_fwhm_arcsec(n_wavelength_obs))
  psf_fwhm_arcsec = 0.0_r2
  if (switch_psf == 1) then
    nf = count(psf_fwhm > 0.0_r2)
    if (nf == 1) then
      psf_fwhm_arcsec(:) = psf_fwhm(1)
    else if (nf == n_wavelength_obs) then
      psf_fwhm_arcsec(:) = psf_fwhm(1:nf)
    else
      call dms_fail("&image psf_fwhm: give one value or one value per wavelength");  return
    end if
  end if
  tel_diameter = telescope_diameter;  mask_radius_arcsec = mask_radius
  select case (lower(image_unit))
  case ("jy/pixel", "jy/pix");     switch_image_unit = 1
  case ("mjy/arcsec2");            switch_image_unit = 2
  case ("jy/arcsec2");             switch_image_unit = 3
  case default;  call dms_fail("&image image_unit must be 'Jy/pixel', 'mJy/arcsec2' or 'Jy/arcsec2'");  return
  end select

  !------------------------------------------------------------------------------------------------
  ! &dust
  !------------------------------------------------------------------------------------------------
  n_grain_components = n_components
  if (n_components < 1 .or. n_components > mc) then
    call dms_fail("&dust n_components must be between 1 and "//integer_to_string(mc));  return
  end if
  allocate(dustfile_name(n_components), particle_bulk_density(n_components), abun(n_components))
  dustfile_name(:)         = optical_data(1:n_components)
  particle_bulk_density(:) = bulk_density(1:n_components)
  abun(:)                  = abundance(1:n_components)

  select case (lower(size_distribution))
  case ("single");               switch_singlegrain_or_Dist = 1
  case ("powerlaw", "power-law"); switch_singlegrain_or_Dist = 2
  case ("file");                 switch_singlegrain_or_Dist = 3
  case default;  call dms_fail("&dust size_distribution must be 'single', 'powerlaw' or 'file'");  return
  end select
  grain_radius_min = a_min;  grain_radius_max = a_max;  q_size = q;  n_grain_radius = n_sizes
  if (switch_singlegrain_or_Dist == 1) then
    grain_radius_max = a_min;  n_grain_radius = 1;  q_size = 0.0_r2
  end if
  sizedist_file = size_file

  switch_dust_mix = merge(1, 0, emt)
  emt_matrix_file = emt_matrix;  emt_inclusion_file = emt_inclusion;  emt_output_file = emt_output
  vol_frac = emt_volume_fraction
  select case (lower(emt_rule))
  case ("maxwell-garnett", "mg");  switch_dust_mix_method = 1
  case ("bruggeman");              switch_dust_mix_method = 2
  case ("lichtenecker");           switch_dust_mix_method = 3
  case ("looyenga");               switch_dust_mix_method = 4
  case default;  call dms_fail("&dust emt_rule unknown: "//trim(emt_rule));  return
  end select

  if (tdep_component > 0) then
    switch_Tdep_optics = 1
    Tdep_comp = tdep_component
    n_Tdep = count(tdep_temperatures > 0.0_r2)
    allocate(Tdep_T(n_Tdep), Tdep_file(n_Tdep))
    Tdep_T(:)    = tdep_temperatures(1:n_Tdep)
    Tdep_file(:) = tdep_files(1:n_Tdep)
  else
    switch_Tdep_optics = 0;  n_Tdep = 0
  end if

  switch_sublimation = merge(1, 0, sublimation)
  allocate(T_subl(n_components), subl_to_comp(n_components))
  T_subl(:)       = t_sublimation(1:n_components)
  subl_to_comp(:) = sublimation_residue(1:n_components)
  if (switch_sublimation == 0) then
    T_subl = -1.0_r2;  subl_to_comp = 0
  end if
  switch_blowout = merge(1, 0, blowout)
  select case (lower(blowout_release))
  case ("circular");                switch_blowout_release = 1
  case ("pericentre", "pericenter"); switch_blowout_release = 2
  case ("apocentre", "apocenter");   switch_blowout_release = 3
  case default
    call dms_fail("&dust blowout_release must be 'circular', 'pericentre' or 'apocentre'");  return
  end select
  blowout_beta_user = blowout_beta
  if (blowout_beta_user > 0.0_r2 .and. switch_blowout_release /= 1) then
    call dms_fail("&dust: give either blowout_beta or blowout_release = 'pericentre'/'apocentre', not both");  return
  end if

  select case (lower(grain_shape))
  case ("sphere", "mie");  switch_grain_shape = 0
  case ("dhs");            switch_grain_shape = 1
  case default;  call dms_fail("&dust grain_shape must be 'sphere' or 'dhs'");  return
  end select
  pm_dhs_fmax = dhs_fmax;  pm_dhs_nf = dhs_nf
  if (switch_grain_shape == 1) then
    ! f_max = 1 is allowed (Min et al. 2005): the Gauss-Legendre nodes lie inside (0, f_max)
    if (.not. (dhs_fmax > 0.0_r2 .and. dhs_fmax <= 1.0_r2)) then
      call dms_fail("&dust dhs_fmax must satisfy 0 < dhs_fmax <= 1");  return
    end if
    if (dhs_nf < 1 .or. dhs_nf > 100) then
      call dms_fail("&dust dhs_nf must be between 1 and 100");  return
    end if
  end if

  !------------------------------------------------------------------------------------------------
  ! &scattering
  !------------------------------------------------------------------------------------------------
  n_angle2 = n_angles
  select case (lower(phase_function))
  case ("mie");                     switch_phase_function = 1
  case ("hg");                      switch_phase_function = 2
  case ("hg2", "two-component-hg"); switch_phase_function = 3
  case default;  call dms_fail("&scattering phase_function must be 'mie', 'hg' or 'hg2'");  return
  end select
  hg_g1 = g1;  hg_g2 = g2;  hg_w = w;  pol_max = p_max
  if (switch_phase_function == 2) hg_w = 1.0_r2

  if (verbose) then
    print *, " "
    print *, "  Input read from ", trim(fname)
    print *, " "
  end if

contains

  subroutine check(group)
    character(len=*), intent(in) :: group
    ! a missing group (end of file) is allowed; anything else is an error in the group
    if (ios /= 0 .and. .not. is_iostat_end(ios)) then
      call dms_fail("error in namelist group &"//group//" of "//trim(fname)//": "//trim(iomsg_txt))
    end if
  end subroutine check

  pure function lower(str) result(res)
    character(len=*), intent(in) :: str
    character(len=len_trim(adjustl(str))) :: res
    integer :: i, ich
    res = trim(adjustl(str))
    do i = 1, len(res)
      ich = iachar(res(i:i))
      if (ich >= 65 .and. ich <= 90) res(i:i) = achar(ich + 32)
    end do
  end function lower

end subroutine read_input_file
