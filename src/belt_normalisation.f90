subroutine belt_normalisation

  !---------------------------------------------------------------------------------------------------
  ! Several belts: the amplitude of every belt is set such that belt k contains the fraction
  ! belt_mass_frac(k) of the (surviving) dust mass:
  !     M_k = int rho_k(x) mfrac(r) dV,     belt_amp(k) = belt_mass_frac(k) / M_k
  ! (mfrac = surviving mass fraction after sublimation / blow-out). The integral is computed in
  ! spherical disk coordinates as in computing_SED. The overall scaling to the total dust mass
  ! is done after the image / SED integration as before.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  integer, parameter :: n_r = 800, n_b = 48, n_p = 360
  integer  :: k, ir, ib, ip
  real(r2) :: ra, rb, rc, b1, b2, bc, psic, dpsi, dV, X, Y, Zd, r_lo, r_hi, mk

  belt_amp(:) = 1.0_r2
  if (n_belts == 1) return

  dpsi = 2.0_r2 * pi / real(n_p, kind=r2)

  do k = 1, n_belts

    r_lo = belt_rin_eff(k)
    r_hi = belt_rout_eff(k)
    if (switch_vertical == 2) r_hi = r_hi / cos(open_angle)

    mk = 0.0_r2
    !$omp parallel do schedule(static) collapse(2) private(ir, ib, ip, ra, rb, rc, b1, b2, bc, psic, dV, X, Y, Zd) &
    !$omp reduction(+:mk)
    do ir = 1, n_r
      do ib = 1, n_b
        ra = r_lo * (r_hi / r_lo)**(real(ir - 1, kind=r2) / real(n_r, kind=r2))
        rb = r_lo * (r_hi / r_lo)**(real(ir,     kind=r2) / real(n_r, kind=r2))
        rc = 0.5_r2 * (ra + rb)
        b1 = -open_angle + 2.0_r2 * open_angle * real(ib - 1, kind=r2) / real(n_b, kind=r2)
        b2 = -open_angle + 2.0_r2 * open_angle * real(ib,     kind=r2) / real(n_b, kind=r2)
        bc = 0.5_r2 * (b1 + b2)
        dV = (rb**3 - ra**3) / 3.0_r2 * (sin(b2) - sin(b1)) * dpsi
        do ip = 1, n_p
          psic = (real(ip, kind=r2) - 0.5_r2) * dpsi
          X  = rc * cos(bc) * sin(psic)
          Y  = rc * cos(bc) * cos(psic)
          Zd = rc * sin(bc)
          mk = mk + density_belt(X, Y, Zd, k) * mfrac_dist(dist_index(rc)) * dV
        end do
      end do
    end do
    !$omp end parallel do

    if (.not. (mk > 0.0_r2)) then
      call dms_fail("belt "//integer_to_string(k)//" contains no (surviving) dust")
      return
    end if
    belt_amp(k) = belt_mass_frac(k) / mk

  end do

  belt_amp(1:n_belts) = belt_amp(1:n_belts) / maxval(belt_amp(1:n_belts))

  if (verbose) then
    if (verbose) print *, " "
    do k = 1, n_belts
      if (verbose) print *, "    belt ", k, ": mass fraction ", real_to_string(100.0_r2 * belt_mass_frac(k)), " %"
    end do
  end if

end subroutine belt_normalisation
