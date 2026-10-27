!____________________________________________________________________________________________________
!
! Python interface of DMS (wrapped with numpy.f2py; see build_extension.py and dms/model.py)
!____________________________________________________________________________________________________
!
! The subroutines below only use plain types so that f2py can wrap them. A run is started from
! an input file (namelist); all results stay in memory until the next run and are copied into
! NumPy arrays by the get_* routines. Errors do not stop the Python interpreter: ierr /= 0 and
! get_error returns the message.
!
! Array order: images are returned as img(ix, iy, il) (Fortran order); dms.model transposes them
! to NumPy order [il, iy, ix].
!____________________________________________________________________________________________________

subroutine run_model(fname, write_files, quiet, ierr)

  use parameters

  implicit none

  character(len=*), intent(in) :: fname
  integer, intent(in)          :: write_files, quiet
  integer, intent(out)         :: ierr

  call dms_reset
  library_mode = .true.
  force_quiet  = (quiet /= 0)

  call read_input_file(trim(fname))
  if (dms_ierr == 0) then
    if (write_files == 0) write_output = .false.
    call dms_run
  end if
  flush(6)

  ierr = dms_ierr

end subroutine run_model


subroutine reset_model

  call dms_reset

end subroutine reset_model


subroutine get_error(msg)

  use parameters

  implicit none

  character(len=512), intent(out) :: msg

  msg = dms_errmsg

end subroutine get_error


subroutine get_dims(nwl, npix, ngr, ncomp, ndist, ntau, has_image, has_obs, complete)

  use parameters

  implicit none

  integer, intent(out) :: nwl, npix, ngr, ncomp, ndist, ntau, has_image, has_obs, complete

  nwl = n_wavelength_obs;  ngr = n_grain_radius;  ncomp = n_grain_components;  ndist = n_dist
  npix = 0;  has_image = 0;  has_obs = 0;  ntau = 0;  complete = 0
  if (allocated(image_total)) then
    npix = number_of_pixels;  has_image = 1
  end if
  if (have_obs) has_obs = 1
  if (allocated(tau_prof_r)) ntau = size(tau_prof_r)
  if (run_complete) complete = 1

end subroutine get_dims


subroutine get_wavelengths(nwl, lam, star_jy)

  use parameters

  implicit none

  integer, intent(in)                    :: nwl
  real(kind=8), dimension(nwl), intent(out) :: lam, star_jy

  lam = wavelength_obs(1:nwl)
  star_jy = stellar_flux_Jy(1:nwl)

end subroutine get_wavelengths


subroutine get_sed(nwl, sed)

  ! sed(:, 1..6) = thermal, scattered, Q, U, polarised, total (dust + star)  [Jy]

  use parameters

  implicit none

  integer, intent(in)                        :: nwl
  real(kind=8), dimension(nwl, 6), intent(out) :: sed

  sed = 0.0d0
  if (.not. allocated(sed_reem)) return
  sed(:, 1) = sed_reem(1:nwl);  sed(:, 2) = sed_sca(1:nwl);  sed(:, 3) = sed_Q(1:nwl)
  sed(:, 4) = sed_U(1:nwl);     sed(:, 5) = sed_pol(1:nwl);  sed(:, 6) = sed_total(1:nwl)

end subroutine get_sed


subroutine get_image(name, npix, nwl, img, ok)

  ! name: reem, sca, q, u, pol, qphi, uphi, poldeg, total  (+ "_obs" for simulated observations)
  ! always in Jy/pixel (poldeg dimensionless)

  use parameters

  implicit none

  character(len=*), intent(in)                     :: name
  integer, intent(in)                              :: npix, nwl
  real(kind=8), dimension(npix, npix, nwl), intent(out) :: img
  integer, intent(out)                             :: ok
  integer :: il

  img = 0.0d0
  ok  = 0
  if (.not. allocated(image_total)) return
  if (number_of_pixels /= npix .or. n_wavelength_obs /= nwl) return

  ok = 1
  do il = 1, nwl
    select case (trim(name))
    case ("reem");        img(:, :, il) = image_reem(il, :, :)
    case ("sca");         img(:, :, il) = image_sca(il, :, :)
    case ("q");           img(:, :, il) = image_Q(il, :, :)
    case ("u");           img(:, :, il) = image_U(il, :, :)
    case ("pol");         img(:, :, il) = image_pol(il, :, :)
    case ("qphi");        img(:, :, il) = image_Qphi(il, :, :)
    case ("uphi");        img(:, :, il) = image_Uphi(il, :, :)
    case ("poldeg");      img(:, :, il) = image_poldeg(il, :, :)
    case ("total");       img(:, :, il) = image_total(il, :, :)
    case default
      if (.not. have_obs) then
        ok = 0;  return
      end if
      select case (trim(name))
      case ("reem_obs");    img(:, :, il) = obs_reem(il, :, :)
      case ("sca_obs");     img(:, :, il) = obs_sca(il, :, :)
      case ("q_obs");       img(:, :, il) = obs_Q(il, :, :)
      case ("u_obs");       img(:, :, il) = obs_U(il, :, :)
      case ("pol_obs");     img(:, :, il) = obs_pol(il, :, :)
      case ("qphi_obs");    img(:, :, il) = obs_Qphi(il, :, :)
      case ("uphi_obs");    img(:, :, il) = obs_Uphi(il, :, :)
      case ("poldeg_obs");  img(:, :, il) = obs_poldeg(il, :, :)
      case ("total_obs");   img(:, :, il) = obs_total(il, :, :)
      case default;         ok = 0;  return
      end select
    end select
  end do

end subroutine get_image


subroutine get_maps(npix, tau, mass)

  ! line-of-sight geometrical optical depth and dust mass per pixel [kg]

  use parameters

  implicit none

  integer, intent(in)                          :: npix
  real(kind=8), dimension(npix, npix), intent(out) :: tau, mass

  tau = 0.0d0;  mass = 0.0d0
  if (.not. allocated(image_tau)) return
  tau  = image_tau(1, :, :)
  mass = total_mass(:, :)

end subroutine get_maps


subroutine get_scalars(vals)

  ! 1 f_d   2 cross section [au^2]   3 tau_mean   4 peak vertical tau   5 pixel [arcsec]
  ! 6 pixel [au]   7 L* [L_sun]   8 T_eff [K]   9 dust mass [M_sun]   10 max. T difference [K]
  ! 11 run complete (1/0)   12 largest radial optical depth in the mid-plane
  ! 13, 14 smallest and largest blow-out threshold beta_crit (0 if blow-out is off)

  use constants
  use parameters

  implicit none

  real(kind=8), dimension(14), intent(out) :: vals

  vals = 0.0d0
  vals(1) = f_dust
  vals(2) = total_cross_section / au**2
  vals(3) = tau_mean
  vals(4) = tau_vertical_peak
  vals(5) = pixel_scale_arcsec
  vals(6) = pixel_size_m / au
  vals(7) = Lstar / 3.828d26
  vals(8) = Tstar
  vals(9) = disk_mass_Msun
  vals(10) = maxTempDiff
  if (run_complete) vals(11) = 1.0d0
  vals(12) = tau_radial_max
  if (switch_blowout == 1 .and. allocated(beta_crit_dist)) then
    vals(13) = beta_crit_min
    vals(14) = beta_crit_max
  end if

end subroutine get_scalars


subroutine get_temperature(ndist, ngr, ncomp, r_au, temp)

  use constants
  use parameters

  implicit none

  integer, intent(in)                                    :: ndist, ngr, ncomp
  real(kind=8), dimension(ndist), intent(out)            :: r_au
  real(kind=8), dimension(ngr, ncomp, ndist), intent(out) :: temp

  r_au = 0.0d0;  temp = 0.0d0
  if (.not. allocated(T_dist)) return
  r_au = dist_grid / au
  temp = T_dist

end subroutine get_temperature


subroutine get_grains(ngr, ncomp, a, edges, number, beta)

  use parameters

  implicit none

  integer, intent(in)                             :: ngr, ncomp
  real(kind=8), dimension(ngr), intent(out)       :: a
  real(kind=8), dimension(ngr + 1), intent(out)   :: edges
  real(kind=8), dimension(ngr, ncomp), intent(out) :: number, beta

  a = rad;  edges = rad_edges
  number = 0.0d0;  beta = 0.0d0
  if (allocated(number_total)) number = number_total
  if (allocated(beta_rp))      beta   = beta_rp

end subroutine get_grains


subroutine get_tau_profile(ntau, prof)

  ! prof(:, 1..4) = distance [au], azimuthal mean, towards pericentre, towards apocentre

  use parameters

  implicit none

  integer, intent(in)                          :: ntau
  real(kind=8), dimension(ntau, 4), intent(out) :: prof

  prof = 0.0d0
  if (.not. allocated(tau_prof_r)) return
  prof(:, 1) = tau_prof_r;    prof(:, 2) = tau_prof_mean
  prof(:, 3) = tau_prof_peri; prof(:, 4) = tau_prof_apo

end subroutine get_tau_profile


subroutine get_phase(nwl, nbin, angle, fsca, fpol)

  ! scattered / polarised flux of the whole disk per 1-degree bin of the scattering angle [Jy]
  ! (fpol > 0: polarisation perpendicular to the scattering plane)

  use parameters

  implicit none

  integer, intent(in)                          :: nwl, nbin
  real(kind=8), dimension(nbin), intent(out)     :: angle
  real(kind=8), dimension(nwl, nbin), intent(out) :: fsca, fpol
  integer :: ib

  do ib = 1, nbin
    angle(ib) = dble(ib) - 0.5d0
  end do
  fsca = 0.0d0;  fpol = 0.0d0
  if (.not. allocated(phase_sca)) return
  fsca = phase_sca;  fpol = phase_pol

end subroutine get_phase


subroutine get_cube(name, npix, nwl, cube, ok)

  ! 3-d cubes (only with &image save_3d_cubes = .true.): flux of every volume element towards the
  ! observer [Jy], cube(ix, iy, iz, il); name = 'reem' or 'sca'

  use parameters

  implicit none

  character(len=*), intent(in)                           :: name
  integer, intent(in)                                    :: npix, nwl
  real(kind=8), dimension(npix, npix, npix, nwl), intent(out) :: cube
  integer, intent(out)                                   :: ok
  integer :: il

  cube = 0.0d0;  ok = 0
  if (.not. allocated(cube_sca)) return
  if (size(cube_sca, 2) /= npix .or. size(cube_sca, 1) /= nwl) return
  ok = 1
  do il = 1, nwl
    select case (trim(name))
    case ("reem");  cube(:, :, :, il) = cube_reem(il, :, :, :)
    case ("sca");   cube(:, :, :, il) = cube_sca(il, :, :, :)
    case default;   ok = 0;  return
    end select
  end do

end subroutine get_cube


subroutine get_cube_angle(npix, ang, ok)

  use parameters

  implicit none

  integer, intent(in)                                :: npix
  real(kind=8), dimension(npix, npix, npix), intent(out) :: ang
  integer, intent(out)                               :: ok

  ang = 0.0d0;  ok = 0
  if (.not. allocated(cube_angle)) return
  if (size(cube_angle, 1) /= npix) return
  ang = cube_angle;  ok = 1

end subroutine get_cube_angle
