!_______________________________________________________________________________________________|
!                                                                                               |
! 5th module : radiation of one (sub-)volume element                                            |
!_______________________________________________________________________________________________|
!===============================================================================================|
! radiation_point adds the contribution of a volume element (centre xd, yd, z in the detector   |
! frame, volume dV) to the accumulators of the calling pixel.                                   |
!                                                                                               |
! All quantities are per unit solid angle towards the observer and are converted to flux        |
! densities by the caller (x mass scaling / d^2, and lambda^2/c to F_nu):                       |
!                                                                                               |
!   thermal emission  : sum_g N_g Cabs_g B_lambda(T_g)              = m * j_reem(:, id)         |
!   scattered light   : sum_g N_g Csca_g S11_g(theta)/int(S11) F*(r) = m/r^2 * sum K_S11          |
!   pol. scattered    : same with S12 (negative = perpendicular to the scattering plane)        |
!                       Stokes Q = P cos(2 phi), U = P sin(2 phi), phi = atan2(yd, xd)          |
!   with F*(r) = pi R*^2 I*(lambda) / r^2 the stellar flux at the grain (built into K_S11).     |
!                                                                                               |
! Scattering angle (observer at +z, star at the origin): cos(theta) = z / r.                    |
! (The previous version used acos(y/r sin i - z/r cos i), which applied the inclination a       |
!  second time to coordinates that are already in the detector frame.)                          |
!                                                                                               |
! The routine only reads global data: it is called concurrently from OpenMP threads.            |
!===============================================================================================|

module radiation_calculation

  use datatype
  use constants
  use parameters
  use functions
  use orbit_populations

  implicit none

contains

  subroutine radiation_point(xd, yd, z, dV, wsym, reem, sca, qq, uu, mass_acc, sig_acc, &
                             fd_acc, hist, ph_sca, ph_pol)

    real(r2), intent(in)                     :: xd, yd, z, dV
    real(r2), intent(in)                     :: wsym      ! weight for global sums (mirror symmetry)
    real(r2), dimension(:), intent(inout)    :: reem, sca, qq, uu
    real(r2), intent(inout)                  :: mass_acc  ! dust mass (unscaled, unweighted)
    real(r2), intent(inout)                  :: sig_acc   ! geometrical cross section
    real(r2), intent(inout)                  :: fd_acc    ! fractional luminosity (weighted)
    real(r2), dimension(:), intent(inout)    :: hist      ! mass per distance bin (weighted)
    real(r2), dimension(:,:), intent(inout)  :: ph_sca    ! (il, 1-deg bin of the scattering angle):
    real(r2), dimension(:,:), intent(inout)  :: ph_pol    ! scattered / polarised flux (weighted)

    real(r2) :: X, Y, Zd, rho, r, m, inv_r2, theta, t, f, rxy2, c2, s2, w
    integer  :: id, ia, ig, ic, is, ib, ip, np, kx, ob, src
    real(r2), dimension(n_wavelength_obs) :: pol, dsca
    real(r2), dimension(max_belts + 1) :: rho_p
    integer,  dimension(max_belts + 1) :: kx_p, ob_p

    call detector_to_disk(xd, yd, z, X, Y, Zd)
    r  = sqrt(xd*xd + yd*yd + z*z)
    id = dist_index(r)
    ! populations at this point: one (the distance entry) without orbit bins; with orbit bins
    ! (eccentric belts, photosputtering on release orbits) one per belt (orbit_populations.f90)
    call point_populations(X, Y, Zd, id, np, rho_p, kx_p, ob_p)
    if (np == 0) return
    if (np == 1 .and. rho_p(1) <= 0.0_r2) return

    if (do_sca) then
      theta = acos(max(-1.0_r2, min(1.0_r2, z / r))) / deg2rad
      t     = theta / d_angle
      ia    = int(t) + 1
      f     = t - real(ia - 1, kind=r2)
      if (ia >= n_angle2) then
        ia = n_angle2 - 1
        f  = 1.0_r2
      end if
      pol(:)  = 0.0_r2
      dsca(:) = 0.0_r2
    end if

    do ip = 1, np
    rho = rho_p(ip)
    if (rho <= 0.0_r2) cycle
    kx = kx_p(ip)
    ob = ob_p(ip)
    m  = rho * dV

    ! m is the mass before sublimation / blow-out; the surviving dust mass is m * mfrac_dist(id).
    ! hist keeps m (the kernels j_reem, sig_dist, fl_dist, npm are per kg of initial mass).
    mass_acc = mass_acc + m * mfrac_dist(kx)
    sig_acc  = sig_acc + m * sig_dist(kx)
    hist(kx) = hist(kx) + m * wsym
    fd_acc   = fd_acc + wsym * m * fl_dist(kx) / (4.0_r2 * pi * r * r)

    !-------------------------------------------------------------------------------------------
    ! thermal re-emission
    !-------------------------------------------------------------------------------------------
    if (do_reem) reem(:) = reem(:) + m * j_reem(:, kx)

    !-------------------------------------------------------------------------------------------
    ! scattering and polarisation
    !-------------------------------------------------------------------------------------------
    if (do_sca) then

      inv_r2 = m / (r * r)

      if (uniform_population) then

        dsca(:) = dsca(:) + inv_r2 * ((1.0_r2 - f) * Ksca_sum(:, ia) + f * Ksca_sum(:, ia+1))
        pol(:)  = pol(:)  + inv_r2 * ((1.0_r2 - f) * Kpol_sum(:, ia) + f * Kpol_sum(:, ia+1))

      else if (ob == 0) then

        do ic = 1, n_grain_components
          do ig = 1, n_grain_radius
            w = npm(ig, ic, id)
            if (w <= 0.0_r2) cycle
            is = sp_dist(ig, ic, id)
            w  = w * inv_r2
            dsca(:) = dsca(:) + w * ((1.0_r2 - f) * KS11(:, ia, ig, is) + f * KS11(:, ia+1, ig, is))
            pol(:) = pol(:) + w * ((1.0_r2 - f) * KS12(:, ia, ig, is) + f * KS12(:, ia+1, ig, is))
          end do
        end do

      else

        ! orbit bin: every initial grain ends as its final component at this distance
        do src = 1, n_grain_components
          do ig = 1, n_grain_radius
            w = npm_init(ig, src)
            if (w <= 0.0_r2) cycle
            if (ps_lost_ob(ig, src, ob)) then
              ic = dest_lost(ig, src, id)
            else
              ic = dest_keep(ig, src, id)
            end if
            if (ic <= 0) cycle
            is = sp_dist(ig, ic, id)
            w  = w * inv_r2
            dsca(:) = dsca(:) + w * ((1.0_r2 - f) * KS11(:, ia, ig, is) + f * KS11(:, ia+1, ig, is))
            pol(:) = pol(:) + w * ((1.0_r2 - f) * KS12(:, ia, ig, is) + f * KS12(:, ia+1, ig, is))
          end do
        end do

      end if

    end if

    end do                ! populations

    if (do_sca) then

      sca(:) = sca(:) + dsca(:)

      ! flux as a function of the scattering angle (1-degree bins; pol < 0 = perpendicular)
      ib = min(int(theta) + 1, 180)
      ph_sca(:, ib) = ph_sca(:, ib) + wsym * dsca(:)
      ph_pol(:, ib) = ph_pol(:, ib) - wsym * pol(:)

      rxy2 = xd*xd + yd*yd
      if (rxy2 > 0.0_r2) then
        c2 = (xd*xd - yd*yd) / rxy2
        s2 = 2.0_r2 * xd * yd / rxy2
        qq(:) = qq(:) + pol(:) * c2
        uu(:) = uu(:) + pol(:) * s2
      end if

    end if

  end subroutine radiation_point

end module radiation_calculation
