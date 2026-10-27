subroutine computing_SED

  !---------------------------------------------------------------------------------------------------
  ! Spectral energy distribution without an image (switch_image_or_SED = 2).
  !
  ! Direct quadrature in spherical disk coordinates (r, latitude beta, azimuth psi):
  !     X = r cos(beta) sin(psi),  Y = r cos(beta) cos(psi),  Z = r sin(beta)
  !     dV = (r2^3 - r1^3)/3 * (sin beta2 - sin beta1) * dpsi        (exact cell volume)
  !   r    : n_r log-spaced cells between the pericentre Rmin(1-e) and the apocentre Rmax(1+e)
  !   beta : n_b cells within +- the half opening angle (the disk is a wedge)
  !   psi  : n_p cells over 2 pi (needed for scattering and for eccentric disks)
  ! Every cell centre is transformed to the detector frame and passed to radiation_point, so the
  ! SED includes the same physics (scattering angles, polarisation, sublimation, eccentricity,
  ! orientation) as the images. Unlike the image calculation no cubes are sampled outside the
  ! disk, which makes this much faster than integrating an image.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions
  use radiation_calculation
  use output_routines

  implicit none

  integer, parameter :: n_r = 400, n_p = 180
  integer  :: n_b
  integer, parameter :: ns = 32
  integer  :: jb, jp
  real(r2) :: thc, bs, ps, dVs
  integer  :: ir, ib, ip, il
  real(r2) :: ra, rb, rc, b1, b2, bc, p1, psic, dpsi, dV, X, Y, Zd, xd, yd, z, conv
  real(r2) :: macc, sacc, fd_tot
  real(r2), allocatable, dimension(:) :: reem, sca, qq, uu, hist_tot, redges
  real(r2), allocatable, dimension(:,:) :: ph_sca, ph_pol
  integer(i8) :: t0, t1, rate

  if (verbose) print *, " "
  if (verbose) print *, " 7) SED (direct quadrature)"
  if (verbose) print *, " "

  allocate(reem(n_wavelength_obs), sca(n_wavelength_obs), qq(n_wavelength_obs), uu(n_wavelength_obs))
  allocate(hist_tot(n_dist), redges(0:n_r))
  reem = 0.0_r2; sca = 0.0_r2; qq = 0.0_r2; uu = 0.0_r2
  hist_tot = 0.0_r2; macc = 0.0_r2; sacc = 0.0_r2; fd_tot = 0.0_r2
  allocate(ph_sca(n_wavelength_obs, n_phase_bins), ph_pol(n_wavelength_obs, n_phase_bins))
  ph_sca = 0.0_r2;  ph_pol = 0.0_r2

  do ir = 0, n_r
    redges(ir) = r_in_eff * (r_max_3d / r_in_eff)**(real(ir, kind=r2) / real(n_r, kind=r2))
  end do
  dpsi = 2.0_r2 * pi / real(n_p, kind=r2)
  n_b  = 16
  if (switch_vertical == 2) n_b = 48          ! resolve the Gaussian vertical profile

  call system_clock(t0, rate)

  !$omp parallel do schedule(dynamic, 1) collapse(2) default(shared) &
  !$omp private(ir, ib, ip, ra, rb, rc, b1, b2, bc, p1, psic, dV, X, Y, Zd, xd, yd, z, thc, jb, jp, bs, ps, dVs) &
  !$omp reduction(+:reem, sca, qq, uu, hist_tot, macc, sacc, fd_tot, ph_sca, ph_pol)
  do ir = 1, n_r
    do ib = 1, n_b
      ra  = redges(ir - 1)
      rb  = redges(ir)
      rc  = 0.5_r2 * (ra + rb)
      b1  = -open_angle + 2.0_r2 * open_angle * real(ib - 1, kind=r2) / real(n_b, kind=r2)
      b2  = -open_angle + 2.0_r2 * open_angle * real(ib,     kind=r2) / real(n_b, kind=r2)
      bc  = 0.5_r2 * (b1 + b2)
      dV  = (rb**3 - ra**3) / 3.0_r2 * (sin(b2) - sin(b1)) * dpsi
      do ip = 1, n_p
        p1 = real(ip - 1, kind=r2) * dpsi
        psic = p1 + 0.5_r2 * dpsi
        X  = rc * cos(bc) * sin(psic)
        Y  = rc * cos(bc) * cos(psic)
        Zd = rc * sin(bc)
        call disk_to_detector(X, Y, Zd, xd, yd, z)
        thc = acos(max(-1.0_r2, min(1.0_r2, z / rc)))
        if (do_sca .and. thc < 3.0_r2 * d_angle * deg2rad + max(dpsi, b2 - b1)) then
          ! close to the forward direction: sub-sample the cell in latitude and azimuth so that
          ! the (strongly varying) phase function is resolved
          do jb = 1, ns
            bs = b1 + (b2 - b1) * (real(jb, kind=r2) - 0.5_r2) / real(ns, kind=r2)
            dVs = (rb**3 - ra**3) / 3.0_r2 * (sin(b1 + (b2 - b1) * real(jb, kind=r2) / real(ns, kind=r2)) &
                  - sin(b1 + (b2 - b1) * real(jb - 1, kind=r2) / real(ns, kind=r2))) * dpsi / real(ns, kind=r2)
            do jp = 1, ns
              ps = p1 + dpsi * (real(jp, kind=r2) - 0.5_r2) / real(ns, kind=r2)
              call disk_to_detector(rc * cos(bs) * sin(ps), rc * cos(bs) * cos(ps), rc * sin(bs), xd, yd, z)
              call radiation_point(xd, yd, z, dVs, 1.0_r2, reem, sca, qq, uu, macc, sacc, fd_tot, hist_tot, &
                                   ph_sca, ph_pol)
            end do
          end do
        else
          call radiation_point(xd, yd, z, dV, 1.0_r2, reem, sca, qq, uu, macc, sacc, fd_tot, hist_tot, &
                               ph_sca, ph_pol)
        end if
      end do
    end do
  end do
  !$omp end parallel do

  call system_clock(t1)
  if (verbose) print *, "    ", n_r * n_b * n_p, " volume elements, time: ", &
           real_to_string(real(t1 - t0, kind=r2) / real(rate, kind=r2)), " s"

  allocate(mass_hist(n_dist))
  mass_hist  = hist_tot
  mass_model = sum(mass_hist * mfrac_dist)
  if (.not. (mass_model > 0.0_r2)) then
    call dms_fail("no dust inside the model space (check the disk parameters / sublimation)")
    return
  end if
  factor_scaling_mass = disk_mass / mass_model
  f_dust = factor_scaling_mass * fd_tot

  allocate(sed_reem(n_wavelength_obs), sed_sca(n_wavelength_obs), sed_Q(n_wavelength_obs), &
           sed_U(n_wavelength_obs), sed_pol(n_wavelength_obs), sed_total(n_wavelength_obs))
  do il = 1, n_wavelength_obs
    conv = factor_scaling_mass / dist_star**2 * lambda_obs_m(il)**2 / c * 1.0e26_r2
    sed_reem(il) = reem(il) * conv
    sed_sca(il)  = sca(il)  * conv
    sed_Q(il)    = qq(il)   * conv
    sed_U(il)    = uu(il)   * conv
    ph_sca(il, :) = ph_sca(il, :) * conv
    ph_pol(il, :) = ph_pol(il, :) * conv
  end do
  allocate(phase_sca(n_wavelength_obs, n_phase_bins), phase_pol(n_wavelength_obs, n_phase_bins))
  phase_sca = ph_sca;  phase_pol = ph_pol
  sed_pol   = sqrt(sed_Q**2 + sed_U**2)
  sed_total = sed_reem + sed_sca + stellar_flux_Jy

  if (write_output) then
    call make_dir(results_dir)
    call write_sed("direct quadrature")
    if (do_sca) call write_phase()
  end if

  deallocate(reem, sca, qq, uu, hist_tot, redges)

end subroutine computing_SED
