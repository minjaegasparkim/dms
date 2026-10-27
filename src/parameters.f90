!____________________________________________________________________________________________________
!
! 3. 3rd module : global variables of DMS
!____________________________________________________________________________________________________
!
! All global quantities live here. Only the set-up routines write to them; the inner integration
! routines (radiation_calculation, volume_splitting) only READ them, which keeps the OpenMP
! parallel image/SED loops thread-safe.
!
! Index conventions (fixed throughout the code):
!   il   : observing wavelength            1 .. n_wavelength_obs
!   ilr  : wavelength of the optical data   1 .. n_wavelength_int   ("ri grid")
!   ig   : grain-size bin                   1 .. n_grain_radius
!   ic   : dust component (material)        1 .. n_grain_components
!   is   : optical species                  1 .. n_species
!          (a component has one species, or several if temperature-dependent optical constants
!           are used; species_comp(is) gives its component, species_T(is) its temperature tag)
!   ia   : scattering angle                 1 .. n_angle2
!   id   : distance-grid point              1 .. n_dist
!
! Coordinate systems
!   detector frame (xd, yd, z) : image plane (xd, yd), z = line of sight towards the observer
!   sky frame      (x , y , z) : detector frame rotated by the position angle about z
!   disk frame     (X , Y , Z) : Z = disk normal, Y = line of nodes (projected major axis)
!____________________________________________________________________________________________________

module parameters

  use datatype
  use constants

  implicit none

  !-----------------------------------------------------------------------------------------------
  ! 0. Directories (relative to the working directory; can be changed in the input file)
  !-----------------------------------------------------------------------------------------------
  character(len=512) :: ri_dir      = "data/optical_constants/"   ! optical constants (n, k)
  character(len=512) :: results_dir = "results/"                  ! output files
  character(len=512) :: miex_dir    = "cache/"                    ! saved optical properties

  !-----------------------------------------------------------------------------------------------
  ! 0b. Run control (used by the command-line program and by the Python interface)
  !-----------------------------------------------------------------------------------------------
  integer            :: dms_ierr = 0                  ! /= 0 after an error
  character(len=512) :: dms_errmsg = ""
  logical            :: library_mode = .false.        ! .true.: errors return instead of stopping
  logical            :: verbose      = .true.         ! screen output
  logical            :: force_quiet  = .false.        ! Python: suppress screen output
  logical            :: write_output = .true.         ! write result files
  logical            :: run_complete = .false.

  !-----------------------------------------------------------------------------------------------
  ! 1. Star
  !-----------------------------------------------------------------------------------------------
  integer            :: switch_BB_or_readingspec          ! 1 black body, 2 spectrum file, 3 model grid
  character(len=512) :: spec_star                         ! spectrum file (relative to the working dir.)
  character(len=512) :: stellar_grid_dir                  ! directory of a model-atmosphere grid
  real(r2)           :: logg_star                         ! log10(g [cm/s^2]) for the model grid
  real(r2)           :: Tstar                             ! effective temperature [K]
  real(r2)           :: Rstar, Rstar_Rsun                 ! stellar radius [m], [R_sun]
  real(r2)           :: Mstar, Mstar_Msun                 ! stellar mass   [kg], [M_sun]
  real(r2)           :: Lstar                             ! stellar luminosity [W] (computed)
  real(r2)           :: dist_star, dist_star_pc           ! distance observer-star [m], [pc]
  integer            :: n_spec                            ! number of points of the spectrum
  real(r2), allocatable, dimension(:) :: spec_lambda      ! [m]
  real(r2), allocatable, dimension(:) :: spec_I           ! stellar surface intensity [W/m^2/m/sr]

  !-----------------------------------------------------------------------------------------------
  ! 2. Disk geometry and density distribution
  !-----------------------------------------------------------------------------------------------
  integer  :: switch_density_fn                            ! 1 power law, 2 Ertel peak, 3 Gaussian,
                                                           ! 4 smooth broken power law
  real(r2) :: alpha, alpha_out                             ! exponents
  real(r2) :: R_start, R_peak, R_c                         ! [m]
  real(r2) :: mu, sigma                                    ! [m]
  real(r2) :: Rmin, Rmax, Rmin_au, Rmax_au                 ! inner/outer semi-major axis [m], [au]
  real(r2) :: open_angle, open_angle_deg                   ! half opening angle
  real(r2) :: disk_mass, disk_mass_Msun                    ! dust mass [kg], [M_sun]
  real(r2) :: inc, inc_deg                                 ! inclination (0 face-on, 90 edge-on)
  real(r2) :: pa, pa_deg                                   ! position angle of the major axis,
                                                           ! counter-clockwise from +yd
  real(r2) :: ecc                                          ! eccentricity of the confocal orbits
  real(r2) :: omega_peri, omega_peri_deg                   ! argument of pericentre (from +Y to +X)

  ! belts (rings): each belt has its own radial law, radii, eccentricity and pericentre;
  ! inclination, position angle and the vertical structure are common to all belts
  integer, parameter :: max_belts = 10
  integer  :: n_belts = 1
  integer,  dimension(max_belts) :: belt_law = 1
  real(r2), dimension(max_belts) :: belt_alpha = 0.0_r2, belt_alpha_out = 0.0_r2, belt_Rc = 0.0_r2, &
                                    belt_Rstart = 0.0_r2, belt_Rpeak = 0.0_r2, belt_mu = 0.0_r2, &
                                    belt_sigma = 0.0_r2                                   ! [m] for radii
  real(r2), dimension(max_belts) :: belt_Rmin = 0.0_r2, belt_Rmax = 0.0_r2               ! [m]
  real(r2), dimension(max_belts) :: belt_ecc = 0.0_r2, belt_omega = 0.0_r2               ! omega [rad]
  real(r2), dimension(max_belts) :: belt_cos_om = 1.0_r2, belt_sin_om = 0.0_r2
  real(r2), dimension(max_belts) :: belt_mass_frac = 1.0_r2                              ! mass fraction
  real(r2), dimension(max_belts) :: belt_amp = 1.0_r2                                    ! normalisation
  real(r2), dimension(max_belts) :: belt_rin_eff = 0.0_r2, belt_rout_eff = 0.0_r2        ! [m]
  real(r2), dimension(max_belts) :: belt_width = 0.0_r2   ! smallest radial width [m] (sampling)

  ! user-defined density (density law 5, see user_density.f90)
  real(r2), dimension(10) :: user_par = 0.0_r2
  logical  :: user_density_symmetric = .false.            ! mirror symmetry yd -> -yd may be used

  ! vertical structure: 1 wedge (constant density within +- open_angle),
  !                     2 Gaussian with scale height h = h_over_r * varpi, cut at gauss_cut * h
  integer  :: switch_vertical = 1
  real(r2) :: h_over_r = 0.05_r2, gauss_cut = 4.0_r2

  ! derived geometry
  real(r2) :: sin_open, tan_open
  real(r2) :: cos_inc, sin_inc, cos_pa, sin_pa, cos_om, sin_om
  real(r2) :: r_in_eff, r_out_eff                          ! smallest pericentre, largest apocentre [m]
  real(r2) :: r_max_3d                                     ! largest distance of dust from the star [m]
  real(r2) :: R_field                                      ! half width of the image / model box [m]

  !-----------------------------------------------------------------------------------------------
  ! 3. Images and dynamic sampling
  !-----------------------------------------------------------------------------------------------
  integer  :: number_of_pixels
  real(r2) :: maxdelT                                      ! max. temperature difference in a cell
  integer  :: Nmax                                         ! max. number of cell splittings
  real(r2) :: pixel_size_m, pixel_scale_arcsec
  logical  :: use_mirror_symmetry                          ! compute half image, mirror in y

  ! simulated observations
  integer  :: switch_psf = 0                               ! 0 none, 1 Gaussian FWHM, 2 diffraction
  real(r2), allocatable, dimension(:) :: psf_fwhm_arcsec   ! (il) FWHM [arcsec]
  real(r2) :: tel_diameter = 0.0_r2                        ! telescope diameter [m] (switch_psf = 2)
  real(r2) :: mask_radius_arcsec = 0.0_r2                  ! coronagraph / inner mask radius [arcsec]
  integer  :: switch_image_unit = 1                        ! 1 Jy/pixel, 2 mJy/arcsec^2, 3 Jy/arcsec^2

  !-----------------------------------------------------------------------------------------------
  ! 4. Wavelengths
  !-----------------------------------------------------------------------------------------------
  integer  :: switch_wavelength_obs                        ! 1 list, 2 log distribution
  integer  :: n_wavelength_obs
  real(r2) :: wavelength_obs_min, wavelength_obs_max       ! [micron]
  real(r2), allocatable, dimension(:) :: wavelength_obs    ! [micron]
  real(r2), allocatable, dimension(:) :: lambda_obs_m      ! [m]

  integer  :: n_wavelength_int                             ! optical-data ("ri") grid
  real(r2), allocatable, dimension(:) :: wavelength_int    ! [micron]
  real(r2), allocatable, dimension(:) :: dlambda_int       ! trapezoidal weights [m]

  !-----------------------------------------------------------------------------------------------
  ! 5. Dust components, optical species, grain sizes
  !-----------------------------------------------------------------------------------------------
  integer :: n_grain_components
  character(len=512), allocatable, dimension(:) :: dustfile_name
  real(r2), allocatable, dimension(:) :: particle_bulk_density   ! [g/cm^3]
  real(r2), allocatable, dimension(:) :: abun                    ! relative abundance (fraction)

  integer :: n_species
  integer,  allocatable, dimension(:) :: species_comp            ! component of a species
  real(r2), allocatable, dimension(:) :: species_T               ! temperature tag [K] (-1: none)
  character(len=512), allocatable, dimension(:) :: species_file
  integer,  allocatable, dimension(:) :: comp_first_species, comp_n_species

  real(r2), allocatable, dimension(:,:) :: n_ri, k_ri            ! (ilr, is) on the ri grid

  integer  :: switch_singlegrain_or_Dist                         ! 1 single, 2 power law, 3 from file
  character(len=512) :: sizedist_file                            ! a [micron], relative number
  integer  :: n_grain_radius
  real(r2) :: grain_radius_min, grain_radius_max                 ! [micron]
  real(r2) :: q_size                                             ! n(a) da ~ a^(-q_size) da
  real(r2), allocatable, dimension(:)   :: rad                   ! representative radius [micron]
  real(r2), allocatable, dimension(:)   :: rad_edges             ! bin edges (0:n) [micron]
  real(r2), allocatable, dimension(:)   :: size_weight           ! number in bin (arbitrary norm.)
  real(r2), allocatable, dimension(:,:) :: grain_mass            ! (ig, ic) [kg]

  !-----------------------------------------------------------------------------------------------
  ! 6. Effective medium theory
  !-----------------------------------------------------------------------------------------------
  integer            :: switch_dust_mix, switch_dust_mix_method
  character(len=512) :: emt_matrix_file, emt_inclusion_file, emt_output_file
  real(r2)           :: vol_frac

  !-----------------------------------------------------------------------------------------------
  ! 7. Temperature-dependent optical constants
  !-----------------------------------------------------------------------------------------------
  integer :: switch_Tdep_optics                                  ! 0 no, 1 yes
  integer :: Tdep_comp                                           ! component that uses them
  integer :: n_Tdep                                              ! number of data sets
  real(r2), allocatable, dimension(:) :: Tdep_T                  ! temperatures of the data sets
  character(len=512), allocatable, dimension(:) :: Tdep_file

  !-----------------------------------------------------------------------------------------------
  ! 8. Optical properties
  !-----------------------------------------------------------------------------------------------
  ! on the ri grid (temperature calculation, radiation pressure)
  real(r2), allocatable, dimension(:,:,:) :: Qabs_ri, Qpr_ri      ! (ilr, ig, is)
  real(r2), allocatable, dimension(:,:)   :: Qstar_mean          ! <Qabs> weighted by the star (ig, is)
  ! at the observing wavelengths
  real(r2), allocatable, dimension(:,:,:)   :: Cabs, Csca         ! (il, ig, is) [m^2]
  real(r2), allocatable, dimension(:,:,:,:) :: S11_norm, S12_norm ! (il, ia, ig, is) S11/(k^2 Csca) [1/sr]
  logical :: optics_from_cache = .false.
  real(r2), allocatable, dimension(:,:,:,:) :: KS11, KS12         ! (il, ia, ig, is):
                                                                 ! Csca * S11/int(S11) * pi R*^2 I*
                                                                 ! [W/m/sr per unit (1/r^2)]
  integer  :: n_angle2, n_angle                                  ! scattering angles (0..180 deg)
  real(r2) :: d_angle                                            ! step [deg]
  real(r2), allocatable, dimension(:) :: x_angle                 ! [deg]

  integer  :: switch_phase_function                              ! 1 Mie, 2 HG, 3 two-component HG
  real(r2) :: hg_g1, hg_g2, hg_w, pol_max

  real(r2), allocatable, dimension(:) :: star_I_obs              ! I*(lambda_obs) [W/m^2/m/sr]
  real(r2), allocatable, dimension(:) :: stellar_flux_Jy         ! stellar flux at observer [Jy]

  !-----------------------------------------------------------------------------------------------
  ! 9. Temperature, sublimation, blow-out, particle numbers
  !-----------------------------------------------------------------------------------------------
  integer  :: n_temp_table = 3000                                ! temperature table 1 K .. 5000 K
  integer  :: n_dist       = 20000                               ! distance grid
  real(r2) :: dist_lo, dist_hi                                   ! [m]
  real(r2) :: dist_step                                          ! step in ln(r)
  real(r2), allocatable, dimension(:)     :: dist_grid           ! [m]
  real(r2), allocatable, dimension(:,:,:) :: T_dist              ! (ig, ic, id) grain temperature [K]
  integer,  allocatable, dimension(:,:,:) :: sp_dist             ! (ig, ic, id) optical species used
  real(r2), allocatable, dimension(:)     :: T1_dist             ! T of grain 1 of component 1

  integer  :: switch_sublimation
  real(r2), allocatable, dimension(:) :: T_subl                  ! (ic) [K]
  integer,  allocatable, dimension(:) :: subl_to_comp            ! (ic) 0: vanishes, n: becomes n

  integer  :: switch_blowout

  ! blow-out threshold: grains with beta >= beta_crit are removed. Parent bodies on an orbit with
  ! eccentricity e release grains that are unbound for beta >= r/(2a) = (1-e^2)/(2(1+e cos f)),
  ! i.e. 1/2 (circular), (1-e)/2 (release at pericentre), (1+e)/2 (release at apocentre).
  !   switch_blowout_release: 1 circular (1/2), 2 pericentre, 3 apocentre (e of the belt)
  !   blowout_beta_user     : > 0 replaces the threshold everywhere
  integer  :: switch_blowout_release = 1
  real(r2) :: blowout_beta_user = -1.0_r2
  real(r2), allocatable, dimension(:) :: beta_crit_dist          ! (id) threshold at every distance
  real(r2) :: beta_crit_min = 0.5_r2, beta_crit_max = 0.5_r2

  ! grain shape: 0 = homogeneous spheres (Mie); 1 = distribution of hollow spheres (DHS, Min et al.
  ! 2005): vacuum volume fraction f uniform in [0, dhs_fmax], dhs_nf Gauss-Legendre points; the
  ! material volume of a grain (and so its mass) does not depend on f
  integer  :: switch_grain_shape
  real(r2) :: dhs_fmax
  integer  :: dhs_nf
  real(r2), allocatable, dimension(:,:) :: beta_rp               ! (ig, ic) F_rad / F_grav

  real(r2), allocatable, dimension(:,:,:) :: npm                 ! (ig, ic, id) particles per kg
  logical  :: uniform_population                                 ! npm independent of distance

  ! kernels per unit dust mass
  real(r2), allocatable, dimension(:,:) :: j_reem                ! (il, id) thermal emission [W/m/sr/kg]
  real(r2), allocatable, dimension(:)   :: fl_dist               ! (id) sum N pi a^2 <Qabs>*  [m^2/kg]
  real(r2), allocatable, dimension(:)   :: sig_dist              ! (id) sum N pi a^2          [m^2/kg]
  real(r2), allocatable, dimension(:)   :: mfrac_dist            ! (id) surviving mass fraction
                                                                 !      (< 1 after sublimation/blow-out)
  real(r2), allocatable, dimension(:,:) :: Ksca_sum, Kpol_sum    ! (il, ia) fast path [W/m/sr m^2/kg]

  !-----------------------------------------------------------------------------------------------
  ! 10. Results
  !-----------------------------------------------------------------------------------------------
  real(r2), allocatable, dimension(:,:,:) :: image_reem, image_sca, image_Q, image_U, &
                                             image_pol, image_total, image_tau   ! (il, ix, iy)
  real(r2), allocatable, dimension(:,:,:) :: image_Qphi, image_Uphi, image_poldeg  ! (il, ix, iy)
  real(r2), allocatable, dimension(:,:,:) :: obs_total, obs_Q, obs_U, obs_pol, obs_Qphi, obs_Uphi, &
                                             obs_poldeg, obs_reem, obs_sca         ! convolved/masked
  logical :: have_obs = .false.
  ! flux as a function of the scattering angle (1-degree bins) and 3-d emission cubes
  integer, parameter :: n_phase_bins = 180
  real(r2), allocatable, dimension(:,:)     :: phase_sca, phase_pol      ! (il, bin) [Jy per bin]
  logical  :: save_cubes = .false.
  real(r2), allocatable, dimension(:,:,:,:) :: cube_reem, cube_sca       ! (il, ix, iy, iz) [Jy]
  real(r2), allocatable, dimension(:,:,:)   :: cube_angle                ! (ix, iy, iz) [deg]
  real(r2), allocatable, dimension(:,:)   :: total_mass                          ! (ix, iy) [kg]
  real(r2), allocatable, dimension(:)     :: tau_prof_r, tau_prof_mean, tau_prof_peri, tau_prof_apo
  real(r2), allocatable, dimension(:)     :: sed_reem, sed_sca, sed_Q, sed_U, sed_pol, sed_total
  real(r2), allocatable, dimension(:)     :: mass_hist                           ! (id) [kg, unscaled]
  real(r2) :: mass_model                                         ! unscaled model mass
  real(r2) :: factor_scaling_mass                                ! disk_mass / mass_model
  real(r2) :: f_dust                                             ! fractional luminosity L_d/L*
  real(r2) :: total_cross_section                                ! sum N pi a^2 [m^2]
  real(r2) :: tau_mean                                           ! cross section / disk area
  real(r2) :: tau_vertical_peak                                  ! max. azimuthally averaged tau_perp
  real(r2) :: tau_radial_max                                     ! max. radial (mid-plane) tau from the star
  real(r2), allocatable, dimension(:,:) :: number_total          ! (ig, ic) number of grains
  real(r2) :: maxTempDiff

  !-----------------------------------------------------------------------------------------------
  ! 11. Switches
  !-----------------------------------------------------------------------------------------------
  integer :: switch_image_or_SED          ! 1 image, 2 SED, 3 both
  integer :: switch_reem_or_scat          ! 1 re-emission, 2 scattering, 3 both, 4 polarisation
  logical :: do_reem, do_sca
  integer :: switch_precal                ! 1 load saved optical properties
  integer :: switch_miex_save
  integer :: switch_temperature_distribution
  integer :: switch_grainsize_distribution
  integer :: switch_file_extensions       ! 0 none, 1 name + parameters, 2 name
  integer :: switch_output_format         ! 1 ASCII, 2 FITS, 3 both
  integer :: switch_add_star              ! add the stellar flux to the central pixel

  character(len=512) :: file_extension_name
  character(len=512) :: file_saved_miex_name

  real(r2) :: var                         ! scratch for input

end module parameters
