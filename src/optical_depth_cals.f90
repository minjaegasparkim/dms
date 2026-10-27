subroutine optical_depth_cals

  !---------------------------------------------------------------------------------------------------
  ! Observational proxies of the dust amount (new)
  !
  !   f_d = L_dust / L*  (fractional luminosity, computed during the integration):
  !         f_d = sum over grains  N pi a^2 <Qabs>* / (4 pi r^2),
  !         <Qabs>* = int Qabs I* dlambda / int I* dlambda
  !
  !   geometrical optical depth  tau = int n pi a^2 ds   (Qext = 1)
  !     * vertical (face-on) optical depth tau_perp(varpi, psi) = int rho(X,Y,Z) sigma/m dZ,
  !       written as azimuthal mean and, for eccentric disks, towards pericentre and apocentre
  !     * total cross section  sum N pi a^2  and the mean optical depth
  !       tau_mean = cross section / disk area,  area = pi (Rmax^2 - Rmin^2) sqrt(1 - e^2)
  !     * the line-of-sight optical depth map image_tau is written by computing_images
  !
  !   number of grains of every size and component, written with the cross sections and beta
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  integer, parameter :: n_w = 400, n_z = 200
  integer  :: n_psi, iw, ipsi, iz, id, ic, ig, u, i_peri, i_apo, k
  real(r2) :: varpi, psi, zmax, dz, Zd, X, Y, r, tau, area, dpsi
  real(r2), allocatable, dimension(:)   :: w_grid, tau_mean_az, tau_peri, tau_apo
  real(r2), allocatable, dimension(:,:) :: tau_map

  if (verbose) print *, " "
  if (verbose) print *, " 8) Fractional luminosity and optical depth"
  if (verbose) print *, " "

  !------------------------------------------------------------------------------------------------
  ! totals
  !------------------------------------------------------------------------------------------------

  total_cross_section = factor_scaling_mass * sum(mass_hist * sig_dist)
  area = 0.0_r2
  do k = 1, n_belts
    area = area + pi * (belt_Rmax(k)**2 - belt_Rmin(k)**2) * sqrt(1.0_r2 - belt_ecc(k)**2)
  end do
  tau_mean = total_cross_section / area

  allocate(number_total(n_grain_radius, n_grain_components))
  do ic = 1, n_grain_components
    do ig = 1, n_grain_radius
      number_total(ig, ic) = factor_scaling_mass * sum(mass_hist(:) * npm(ig, ic, :))
    end do
  end do

  !------------------------------------------------------------------------------------------------
  ! vertical optical depth
  !------------------------------------------------------------------------------------------------

  n_psi = 1
  if (any(belt_ecc(1:n_belts) > 0.0_r2) .or. any(belt_law(1:n_belts) == 5)) n_psi = 72
  dpsi = 2.0_r2 * pi / real(n_psi, kind=r2)

  allocate(w_grid(n_w), tau_map(n_w, n_psi), tau_mean_az(n_w), tau_peri(n_w), tau_apo(n_w))

  !$omp parallel do schedule(dynamic) private(iw, ipsi, iz, varpi, psi, zmax, dz, Zd, X, Y, r, id, tau)
  do iw = 1, n_w
    varpi = r_in_eff * (r_out_eff / r_in_eff)**((real(iw, kind=r2) - 0.5_r2) / real(n_w, kind=r2))   ! log grid
    w_grid(iw) = varpi
    do ipsi = 1, n_psi
      psi  = omega_peri + real(ipsi - 1, kind=r2) * dpsi          ! ipsi = 1: pericentre direction
      X    = varpi * sin(psi)
      Y    = varpi * cos(psi)
      zmax = tan_open * varpi
      dz   = 2.0_r2 * zmax / real(n_z, kind=r2)
      tau  = 0.0_r2
      do iz = 1, n_z
        Zd = -zmax + (real(iz, kind=r2) - 0.5_r2) * dz
        r  = sqrt(varpi * varpi + Zd * Zd)
        id = dist_index(r)
        tau = tau + density_disk(X, Y, Zd) * sig_dist(id) * dz
      end do
      tau_map(iw, ipsi) = tau * factor_scaling_mass
    end do
  end do
  !$omp end parallel do

  !------------------------------------------------------------------------------------------------
  ! radial geometrical optical depth in the mid-plane, seen from the star (optically thin check)
  !------------------------------------------------------------------------------------------------
  tau_radial_max = 0.0_r2
  !$omp parallel do schedule(dynamic) private(ipsi, psi, iz, varpi, dz, X, Y, tau) reduction(max:tau_radial_max)
  do ipsi = 1, 72
    psi = real(ipsi - 1, kind=r2) * 2.0_r2 * pi / 72.0_r2
    tau = 0.0_r2
    do iz = 1, 2000
      varpi = r_in_eff * (r_max_3d / r_in_eff)**((real(iz, kind=r2) - 0.5_r2) / 2000.0_r2)
      dz    = varpi * log(r_max_3d / r_in_eff) / 2000.0_r2
      X = varpi * sin(psi);  Y = varpi * cos(psi)
      tau = tau + density_disk(X, Y, 0.0_r2) * sig_dist(dist_index(varpi)) * dz
    end do
    tau_radial_max = max(tau_radial_max, tau * factor_scaling_mass)
  end do
  !$omp end parallel do

  if (verbose .and. tau_radial_max > 0.1_r2) then
    print *, " "
    print *, "    WARNING: the radial optical depth in the mid-plane reaches ", &
             real_to_scientific_string(tau_radial_max), "."
    print *, "    DMS assumes an optically thin disk (no extinction of the stellar light, single"
    print *, "    scattering). Temperatures, fluxes and f_d are not reliable for this dust mass."
  end if

  i_peri = 1
  i_apo  = n_psi / 2 + 1
  do iw = 1, n_w
    tau_mean_az(iw) = sum(tau_map(iw, :)) / real(n_psi, kind=r2)
    tau_peri(iw)    = tau_map(iw, i_peri)
    tau_apo(iw)     = tau_map(iw, i_apo)
  end do
  tau_vertical_peak = maxval(tau_mean_az)

  allocate(tau_prof_r(n_w), tau_prof_mean(n_w), tau_prof_peri(n_w), tau_prof_apo(n_w))
  tau_prof_r = w_grid / au;  tau_prof_mean = tau_mean_az;  tau_prof_peri = tau_peri;  tau_prof_apo = tau_apo

  if (.not. write_output) then
    deallocate(w_grid, tau_map, tau_mean_az, tau_peri, tau_apo)
    return
  end if

  call make_dir(results_dir)
  open(newunit=u, file=trim(results_dir)//"optical_depth"//file_suffix(0)//".dat", action="write", &
       status="replace", form="formatted")
  write(u, '(a)') "# vertical (face-on) geometrical optical depth tau_perp = int n pi a^2 dZ"
  write(u, '(a, es15.6e3)') "# fractional luminosity f_d = L_dust/L* = ", f_dust
  write(u, '(a, es15.6e3, a)') "# total geometrical cross section = ", total_cross_section / au**2, " au^2"
  write(u, '(a, es15.6e3)') "# mean optical depth (cross section / disk area) = ", tau_mean
  write(u, '(a, es15.6e3)') "# largest radial optical depth in the mid-plane (from the star) = ", tau_radial_max
  write(u, '(a)') "# 1 distance from the star in the disk plane [au]  2 azimuthal mean" // &
                  "  3 towards pericentre  4 towards apocentre"
  do iw = 1, n_w
    write(u, '(4es15.6e3)') w_grid(iw) / au, tau_mean_az(iw), tau_peri(iw), tau_apo(iw)
  end do
  close(u)

  !------------------------------------------------------------------------------------------------
  ! grain size distribution
  !------------------------------------------------------------------------------------------------

  if (switch_grainsize_distribution == 1) then
    open(newunit=u, file=trim(results_dir)//"grain_size_distribution"//file_suffix(0)//".dat", &
         action="write", status="replace", form="formatted")
    write(u, '(a)') "# number of grains in the whole disk (after sublimation / blow-out)"
    write(u, '(a)') "# 1 component  2 a [micron]  3 a_low  4 a_high  5 N  6 N pi a^2 [au^2]" // &
                    "  7 mass [kg]  8 beta"
    do ic = 1, n_grain_components
      do ig = 1, n_grain_radius
        write(u, '(i4, 7es15.6e3)') ic, rad(ig), rad_edges(ig-1), rad_edges(ig), number_total(ig, ic), &
              number_total(ig, ic) * pi * (rad(ig) * micron)**2 / au**2, &
              number_total(ig, ic) * grain_mass(ig, ic), beta_rp(ig, ic)
      end do
    end do
    close(u)
  end if

  deallocate(w_grid, tau_map, tau_mean_az, tau_peri, tau_apo)

end subroutine optical_depth_cals
