subroutine observation_products

  !---------------------------------------------------------------------------------------------------
  ! 1. Simulated observations (only if a PSF and/or an inner mask is requested):
  !    the images (thermal, scattered, Stokes Q, U, total incl. the star if added) are convolved
  !    with a Gaussian PSF (FWHM given, or diffraction limit FWHM = 1.029 lambda/D); polarised
  !    intensity, Q_phi, U_phi and the degree of polarisation are derived from the convolved Q, U;
  !    pixels inside the mask radius are set to NaN.  ->  obs_* arrays, files *_obs
  ! 2. Output of all images (unit chosen in the input), radial profiles and (mode 3) the SED.
  !
  ! Internally all images are kept in Jy/pixel (their sum is the SED).
  !---------------------------------------------------------------------------------------------------

  use, intrinsic :: ieee_arithmetic
  use datatype
  use constants
  use parameters
  use functions
  use output_routines
  use obs_tools

  implicit none

  integer  :: N, ic0, il, ix, iy
  real(r2) :: ufac, fwhm, rmask_pix, nan
  real(r2), allocatable :: sigma_pix(:)
  character(len=:), allocatable :: ustr

  N   = number_of_pixels
  ic0 = (N + 1) / 2

  select case (switch_image_unit)
  case (2)
    ufac = 1.0e3_r2 / pixel_scale_arcsec**2;  ustr = "mJy/arcsec2"
  case (3)
    ufac = 1.0_r2 / pixel_scale_arcsec**2;    ustr = "Jy/arcsec2"
  case default
    ufac = 1.0_r2;                            ustr = "Jy/pixel"
  end select

  !------------------------------------------------------------------------------------------------
  ! 1. simulated observations
  !------------------------------------------------------------------------------------------------

  have_obs = (switch_psf > 0 .or. mask_radius_arcsec > 0.0_r2)

  if (have_obs) then

    allocate(sigma_pix(n_wavelength_obs))
    do il = 1, n_wavelength_obs
      select case (switch_psf)
      case (1)
        fwhm = psf_fwhm_arcsec(il)
      case (2)
        fwhm = 1.029_r2 * lambda_obs_m(il) / tel_diameter * (180.0_r2 / pi * 3600.0_r2)
        psf_fwhm_arcsec(il) = fwhm
      case default
        fwhm = 0.0_r2
      end select
      sigma_pix(il) = fwhm / (2.0_r2 * sqrt(2.0_r2 * log(2.0_r2))) / pixel_scale_arcsec
    end do

    allocate(obs_reem(n_wavelength_obs, N, N), obs_sca(n_wavelength_obs, N, N), &
             obs_Q(n_wavelength_obs, N, N), obs_U(n_wavelength_obs, N, N), obs_total(n_wavelength_obs, N, N), &
             obs_pol(n_wavelength_obs, N, N), obs_Qphi(n_wavelength_obs, N, N), &
             obs_Uphi(n_wavelength_obs, N, N), obs_poldeg(n_wavelength_obs, N, N))
    obs_reem = image_reem;  obs_sca = image_sca;  obs_Q = image_Q;  obs_U = image_U;  obs_total = image_total

    if (switch_psf > 0) then
      if (do_reem) call gauss_convolve(obs_reem, sigma_pix)
      if (do_sca) then
        call gauss_convolve(obs_sca, sigma_pix)
        call gauss_convolve(obs_Q,   sigma_pix)
        call gauss_convolve(obs_U,   sigma_pix)
      end if
      call gauss_convolve(obs_total, sigma_pix)
    end if

    obs_pol = sqrt(obs_Q**2 + obs_U**2)
    call azimuthal_stokes(obs_Q, obs_U, obs_Qphi, obs_Uphi, N)
    call pol_degree(obs_pol, obs_reem + obs_sca, obs_poldeg)

    if (mask_radius_arcsec > 0.0_r2) then
      nan = ieee_value(nan, ieee_quiet_nan)
      rmask_pix = mask_radius_arcsec / pixel_scale_arcsec
      do iy = 1, N
        do ix = 1, N
          if (real((ix - ic0)**2 + (iy - ic0)**2, kind=r2) <= rmask_pix**2) then
            obs_reem(:, ix, iy) = nan;  obs_sca(:, ix, iy) = nan;  obs_total(:, ix, iy) = nan
            obs_Q(:, ix, iy) = nan;     obs_U(:, ix, iy) = nan;    obs_pol(:, ix, iy) = nan
            obs_Qphi(:, ix, iy) = nan;  obs_Uphi(:, ix, iy) = nan; obs_poldeg(:, ix, iy) = nan
          end if
        end do
      end do
    end if

    if (verbose) then
      print *, " "
      if (switch_psf > 0) print *, "    simulated observation: Gaussian PSF, FWHM = ", &
                                   real_to_scientific_string(psf_fwhm_arcsec(1)), " arcsec at ", &
                                   real_to_string(wavelength_obs(1)), " micron"
      if (mask_radius_arcsec > 0.0_r2) print *, "    inner mask radius ", real_to_string(mask_radius_arcsec), " arcsec"
    end if

    deallocate(sigma_pix)

  end if

  !------------------------------------------------------------------------------------------------
  ! 2. output
  !------------------------------------------------------------------------------------------------

  if (.not. write_output) return

  call make_dir(results_dir)

  if (do_reem) call write_image("image_reem", image_reem * ufac, ustr, "thermal re-emission", .true.)
  if (do_sca) then
    call write_image("image_sca",    image_sca  * ufac, ustr, "scattered stellar light", .true.)
    call write_image("image_Q",      image_Q    * ufac, ustr, "Stokes Q (+: along the image x axis)", .true.)
    call write_image("image_U",      image_U    * ufac, ustr, "Stokes U", .true.)
    call write_image("image_pol",    image_pol  * ufac, ustr, "polarised intensity sqrt(Q^2+U^2)", .true.)
    call write_image("image_Qphi",   image_Qphi * ufac, ustr, "azimuthal Stokes Q_phi", .true.)
    call write_image("image_Uphi",   image_Uphi * ufac, ustr, "azimuthal Stokes U_phi", .true.)
    call write_image("image_poldeg", image_poldeg,      "1",  "degree of linear polarisation PI/I", .true.)
  end if
  call write_image("image_total", image_total * ufac, ustr, "total intensity (dust + star if added)", .true.)
  call write_image("image_tau",   image_tau,          "1",  "line-of-sight geometrical optical depth", .false.)

  if (have_obs) then
    if (do_reem) call write_image("image_reem_obs", obs_reem * ufac, ustr, "thermal re-emission, observed", .true.)
    if (do_sca) then
      call write_image("image_sca_obs",    obs_sca  * ufac, ustr, "scattered light, observed", .true.)
      call write_image("image_Q_obs",      obs_Q    * ufac, ustr, "Stokes Q, observed", .true.)
      call write_image("image_U_obs",      obs_U    * ufac, ustr, "Stokes U, observed", .true.)
      call write_image("image_pol_obs",    obs_pol  * ufac, ustr, "polarised intensity, observed", .true.)
      call write_image("image_Qphi_obs",   obs_Qphi * ufac, ustr, "Q_phi, observed", .true.)
      call write_image("image_Uphi_obs",   obs_Uphi * ufac, ustr, "U_phi, observed", .true.)
      call write_image("image_poldeg_obs", obs_poldeg,      "1",  "degree of polarisation, observed", .true.)
    end if
    call write_image("image_total_obs", obs_total * ufac, ustr, "total intensity, observed", .true.)
  end if

  call write_profiles(N, ufac, ustr)

  if (switch_image_or_SED == 3) call write_sed("integrated images")
  if (do_sca) call write_phase()
  if (save_cubes) call write_cubes()

end subroutine observation_products
