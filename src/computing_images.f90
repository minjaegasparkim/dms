subroutine computing_images

  !---------------------------------------------------------------------------------------------------
  ! Images (thermal emission, scattered light, Stokes Q/U, polarised intensity, total) and maps of
  ! the line-of-sight geometrical optical depth and of the dust mass.
  !
  ! The model space is a cube of N x N x N volume elements (N = number of pixels, odd) with
  ! half width R_field = largest distance of dust from the star (apocentre of the outermost orbit). Each pixel column is
  ! integrated along the line of sight; each volume element is sub-sampled by SubVol.
  !
  ! Speed-ups with respect to the previous version
  !   * OpenMP over pixels (dynamic schedule; the work per pixel varies strongly)
  !   * mirror symmetry yd -> -yd when the model allows it (PA = 0 or 180 deg, and e = 0 or the
  !     pericentre on the projected minor axis): only half of the image is computed
  !   * empty volume elements are rejected by exact distance and wedge tests in SubVol
  !   * all optical data are pre-computed kernels (no Mie/temperature work inside the loop)
  !
  ! Normalisation: the density law is unscaled; after the integration the model is scaled so that
  ! the surviving dust mass (after sublimation / blow-out) equals the requested dust mass.
  ! Units: images in Jy/pixel, tau dimensionless, mass map in kg/pixel.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions
  use volume_splitting
  use output_routines
  use obs_tools

  implicit none

  integer  :: N, ic0, ix, iy, iz, ip, npix, iy_first, iym, il, k
  integer, allocatable, dimension(:) :: list_x, list_y
  real(r2) :: xd, yd, cz, h0, wsym, macc, sacc, fdacc, fd_tot, rproj, conv
  real(r2), allocatable, dimension(:)   :: reem, sca, qq, uu, hist_tot, reem0, sca0
  real(r2), allocatable, dimension(:,:) :: ph_sca, ph_pol
  real(r2) :: rc, mem_gb
  real(r2), allocatable, dimension(:,:) :: sig_map
  integer(i8) :: t0, t1, rate

  if (verbose) print *, " "
  if (verbose) print *, " 7) Images"
  if (verbose) print *, " "

  N   = number_of_pixels
  ic0 = (N + 1) / 2

  pixel_size_m       = 2.0_r2 * R_field / real(N, kind=r2)
  pixel_scale_arcsec = (pixel_size_m / au) / dist_star_pc
  h0                 = 0.5_r2 * pixel_size_m

  use_mirror_symmetry = (sin_pa == 0.0_r2)
  do k = 1, n_belts
    if (belt_ecc(k) > 0.0_r2 .and. belt_cos_om(k) /= 0.0_r2) use_mirror_symmetry = .false.
    if (belt_law(k) == 5 .and. .not. user_density_symmetric)  use_mirror_symmetry = .false.
  end do

  allocate(image_reem(n_wavelength_obs, N, N), image_sca(n_wavelength_obs, N, N), &
           image_Q(n_wavelength_obs, N, N),    image_U(n_wavelength_obs, N, N))
  allocate(total_mass(N, N), sig_map(N, N), hist_tot(n_dist))
  image_reem = 0.0_r2; image_sca = 0.0_r2; image_Q = 0.0_r2; image_U = 0.0_r2
  total_mass = 0.0_r2; sig_map = 0.0_r2;   hist_tot = 0.0_r2
  fd_tot = 0.0_r2

  ! list of pixels to compute (pixel columns that cannot touch the disk are skipped)
  iy_first = 1
  if (use_mirror_symmetry) iy_first = ic0
  allocate(list_x(N * N), list_y(N * N))
  npix = 0
  do iy = iy_first, N
    do ix = 1, N
      xd = real(ix - ic0, kind=r2) * pixel_size_m
      yd = real(iy - ic0, kind=r2) * pixel_size_m
      rproj = sqrt(max(abs(xd) - h0, 0.0_r2)**2 + max(abs(yd) - h0, 0.0_r2)**2)
      if (rproj > r_max_3d) cycle
      npix = npix + 1
      list_x(npix) = ix
      list_y(npix) = iy
    end do
  end do

  if (verbose) print *, "    ", N, " x ", N, " pixels, pixel = ", real_to_string(pixel_size_m / au), " au = ", &
           real_to_scientific_string(pixel_scale_arcsec), " arcsec"
  if (use_mirror_symmetry .and. verbose) print *, "    mirror symmetry used (half of the image is computed)"
  if (verbose) print *, "    pixel columns to integrate: ", npix

  allocate(ph_sca(n_wavelength_obs, n_phase_bins), ph_pol(n_wavelength_obs, n_phase_bins))
  ph_sca = 0.0_r2;  ph_pol = 0.0_r2

  ! optional 3-d cubes: flux emitted by every volume element towards the observer
  if (save_cubes) then
    mem_gb = 2.0_r2 * 8.0_r2 * real(n_wavelength_obs, kind=r2) * real(N, kind=r2)**3 / 1.0e9_r2
    if (mem_gb > 8.0_r2) then
      call dms_fail("3-d cubes would need "//real_to_string(mem_gb)// &
                    " GB of memory; use fewer pixels or wavelengths");  return
    end if
    allocate(cube_reem(n_wavelength_obs, N, N, N), cube_sca(n_wavelength_obs, N, N, N), cube_angle(N, N, N))
    cube_reem = 0.0_r2;  cube_sca = 0.0_r2
    do iz = 1, N
      do iy = 1, N
        do ix = 1, N
          rc = sqrt(real((ix-ic0)**2 + (iy-ic0)**2 + (iz-ic0)**2, kind=r2))
          if (rc > 0.0_r2) then
            cube_angle(ix, iy, iz) = acos(real(iz - ic0, kind=r2) / rc) / deg2rad
          else
            cube_angle(ix, iy, iz) = 0.0_r2
          end if
        end do
      end do
    end do
  end if

  call system_clock(t0, rate)

  !$omp parallel default(shared) private(ip, ix, iy, iz, xd, yd, cz, wsym, reem, sca, qq, uu, &
  !$omp                                  macc, sacc, fdacc, reem0, sca0) reduction(+:hist_tot, fd_tot, ph_sca, ph_pol)
  allocate(reem(n_wavelength_obs), sca(n_wavelength_obs), qq(n_wavelength_obs), uu(n_wavelength_obs))
  allocate(reem0(n_wavelength_obs), sca0(n_wavelength_obs))

  !$omp do schedule(dynamic, 1)
  do ip = 1, npix

    ix = list_x(ip)
    iy = list_y(ip)
    xd = real(ix - ic0, kind=r2) * pixel_size_m
    yd = real(iy - ic0, kind=r2) * pixel_size_m

    wsym = 1.0_r2
    if (use_mirror_symmetry .and. iy > ic0) wsym = 2.0_r2

    reem = 0.0_r2; sca = 0.0_r2; qq = 0.0_r2; uu = 0.0_r2
    macc = 0.0_r2; sacc = 0.0_r2; fdacc = 0.0_r2

    do iz = 1, N
      cz = real(iz - ic0, kind=r2) * pixel_size_m
      if (save_cubes) then
        reem0 = reem;  sca0 = sca
      end if
      call SubVol(xd, yd, cz, h0, wsym, reem, sca, qq, uu, macc, sacc, fdacc, hist_tot, ph_sca, ph_pol)
      if (save_cubes) then
        cube_reem(:, ix, iy, iz) = reem - reem0
        cube_sca(:, ix, iy, iz)  = sca  - sca0
      end if
    end do

    image_reem(:, ix, iy) = reem
    image_sca(:, ix, iy)  = sca
    image_Q(:, ix, iy)    = qq
    image_U(:, ix, iy)    = uu
    total_mass(ix, iy)    = macc
    sig_map(ix, iy)       = sacc
    fd_tot                = fd_tot + fdacc

  end do
  !$omp end do

  deallocate(reem, sca, qq, uu, reem0, sca0)
  !$omp end parallel

  call system_clock(t1)
  if (verbose) print *, "    integration time: ", real_to_string(real(t1 - t0, kind=r2) / real(rate, kind=r2)), " s"

  ! mirror (yd -> -yd: Stokes U changes sign)
  if (use_mirror_symmetry) then
    do iy = ic0 + 1, N
      iym = N + 1 - iy
      image_reem(:, :, iym) =  image_reem(:, :, iy)
      image_sca(:, :, iym)  =  image_sca(:, :, iy)
      image_Q(:, :, iym)    =  image_Q(:, :, iy)
      image_U(:, :, iym)    = -image_U(:, :, iy)
      if (save_cubes) then
        cube_reem(:, :, iym, :) = cube_reem(:, :, iy, :)
        cube_sca(:, :, iym, :)  = cube_sca(:, :, iy, :)
      end if
      total_mass(:, iym)    =  total_mass(:, iy)
      sig_map(:, iym)       =  sig_map(:, iy)
    end do
  end if

  !------------------------------------------------------------------------------------------------
  ! normalisation to the dust mass
  !------------------------------------------------------------------------------------------------

  allocate(mass_hist(n_dist))
  mass_hist  = hist_tot
  mass_model = sum(mass_hist * mfrac_dist)
  if (.not. (mass_model > 0.0_r2)) then
    call dms_fail("no dust inside the model space (check the disk parameters / sublimation)")
    return
  end if
  factor_scaling_mass = disk_mass / mass_model
  f_dust = factor_scaling_mass * fd_tot

  if (verbose) print *, "    mass check: sum of the mass map / dust mass = ", &
           real_to_string(factor_scaling_mass * sum(total_mass) / disk_mass)

  do il = 1, n_wavelength_obs
    conv = factor_scaling_mass / dist_star**2 * lambda_obs_m(il)**2 / c * 1.0e26_r2   ! -> Jy
    image_reem(il, :, :) = image_reem(il, :, :) * conv
    image_sca(il, :, :)  = image_sca(il, :, :)  * conv
    image_Q(il, :, :)    = image_Q(il, :, :)    * conv
    image_U(il, :, :)    = image_U(il, :, :)    * conv
    ph_sca(il, :)        = ph_sca(il, :)        * conv
    ph_pol(il, :)        = ph_pol(il, :)        * conv
    if (save_cubes) then
      cube_reem(il, :, :, :) = cube_reem(il, :, :, :) * conv
      cube_sca(il, :, :, :)  = cube_sca(il, :, :, :)  * conv
    end if
  end do
  total_mass = total_mass * factor_scaling_mass

  allocate(image_pol(n_wavelength_obs, N, N), image_total(n_wavelength_obs, N, N), &
           image_tau(1, N, N))
  image_pol   = sqrt(image_Q**2 + image_U**2)
  image_total = image_reem + image_sca
  if (switch_add_star == 1) image_total(:, ic0, ic0) = image_total(:, ic0, ic0) + stellar_flux_Jy(:)
  image_tau(1, :, :) = factor_scaling_mass * sig_map(:, :) / pixel_size_m**2

  ! azimuthal Stokes parameters and degree of polarisation
  allocate(image_Qphi(n_wavelength_obs, N, N), image_Uphi(n_wavelength_obs, N, N), &
           image_poldeg(n_wavelength_obs, N, N))
  call azimuthal_stokes(image_Q, image_U, image_Qphi, image_Uphi, N)
  call pol_degree(image_pol, image_reem + image_sca, image_poldeg)

  !------------------------------------------------------------------------------------------------
  ! SED from the images (spatially integrated)
  !------------------------------------------------------------------------------------------------

  allocate(sed_reem(n_wavelength_obs), sed_sca(n_wavelength_obs), sed_Q(n_wavelength_obs), &
           sed_U(n_wavelength_obs), sed_pol(n_wavelength_obs), sed_total(n_wavelength_obs))
  do il = 1, n_wavelength_obs
    sed_reem(il) = sum(image_reem(il, :, :))
    sed_sca(il)  = sum(image_sca(il, :, :))
    sed_Q(il)    = sum(image_Q(il, :, :))
    sed_U(il)    = sum(image_U(il, :, :))
  end do
  sed_pol   = sqrt(sed_Q**2 + sed_U**2)
  sed_total = sed_reem + sed_sca + stellar_flux_Jy

  allocate(phase_sca(n_wavelength_obs, n_phase_bins), phase_pol(n_wavelength_obs, n_phase_bins))
  phase_sca = ph_sca;  phase_pol = ph_pol

  deallocate(list_x, list_y, sig_map, hist_tot, ph_sca, ph_pol)

end subroutine computing_images
