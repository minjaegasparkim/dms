subroutine Qabs_Qsca_S11_S12_cals

  !========================================================================================|
  ! Optical properties at the observing wavelengths                                        |
  !                                                                                        |
  ! 1. Mie calculation at the observing wavelengths with the refractive index interpolated |
  !    to lambda_obs (instead of a linear interpolation of Qabs, Qsca, S11 between the     |
  !    neighbouring grid wavelengths, as before). Gives Cabs, Csca and the normalised      |
  !    scattering matrix elements                                                          |
  !        S11_norm = S11 / (k^2 Csca),  S12_norm = S12 / (k^2 Csca)   [1/sr]              |
  !    with int S11_norm dOmega = 1 (k = 2 pi / lambda; Bohren & Huffman 1983).             |
  ! 2. Optional replacement of the Mie phase function by a (two-component)                 |
  !    Henyey-Greenstein function                                                          |
  !        p(theta) = w HG(g1, theta) + (1 - w) HG(g2, theta)                              |
  !    and a Rayleigh-like degree of polarisation                                          |
  !        P(theta) = -p_max sin^2(theta) / (1 + cos^2(theta)),  S12_norm = P * p(theta)   |
  !    (Engler et al. 2017, A&A 607, A90). Cabs and Csca stay those of Mie theory.         |
  ! 3. Scattering kernels  K_S11 = Csca * S11_norm * pi R*^2 I*(lambda) used by            |
  !    radiation_point (the stellar flux at distance r is pi R*^2 I* / r^2).               |
  !========================================================================================|

  use datatype
  use constants
  use parameters
  use functions
  use mie_routines

  implicit none

  integer     :: ntot, it, il, ig, is, ia, ier
  real(r2)    :: sizepar, qext, qsca, qabs, qbk, qpr, alb, gsc, nn, kk, norm, ct, ph, pp
  complex(r2) :: ri
  integer     :: mfine, kfwd, nab
  complex(r2), allocatable, dimension(:,:) :: ab
  complex(r2), allocatable, dimension(:) :: sf1, sf2
  real(r2),    allocatable, dimension(:) :: f11, f12
  real(r2), allocatable, dimension(:,:,:) :: Qabs_o, Qsca_o

  if (verbose) print *, " "
  if (verbose) print *, " 4) Optical properties at the observing wavelengths"
  if (verbose) print *, " "

  !========================================================================================
  ! 1. Mie theory at lambda_obs
  !========================================================================================

  if (.not. optics_from_cache) then

    allocate(Cabs(n_wavelength_obs, n_grain_radius, n_species), &
             Csca(n_wavelength_obs, n_grain_radius, n_species))
    allocate(Qabs_o(n_wavelength_obs, n_grain_radius, n_species), &
             Qsca_o(n_wavelength_obs, n_grain_radius, n_species))
    if (do_sca) then
      allocate(S11_norm(n_wavelength_obs, n_angle2, n_grain_radius, n_species), &
               S12_norm(n_wavelength_obs, n_angle2, n_grain_radius, n_species))
    end if

    ntot = n_wavelength_obs * n_grain_radius * n_species

    !$omp parallel do schedule(dynamic, 1) default(shared) &
    !$omp private(it, il, ig, is, nn, kk, ri, sizepar, qext, qsca, qabs, qbk, qpr, alb, gsc, ier, &
    !$omp         norm, mfine, kfwd, nab, ab, sf1, sf2, f11, f12)
    do it = 1, ntot

      il = mod(it - 1, n_wavelength_obs) + 1
      ig = mod((it - 1) / n_wavelength_obs, n_grain_radius) + 1
      is = (it - 1) / (n_wavelength_obs * n_grain_radius) + 1

      nn = interp_lin(log(wavelength_int), n_ri(:, is), log(wavelength_obs(il)))
      kk = interp_lin(log(wavelength_int), k_ri(:, is), log(wavelength_obs(il)))
      ri = cmplx(nn, max(kk, 0.0_r2), kind=r2)

      sizepar = 2.0_r2 * pi * rad(ig) / wavelength_obs(il)

      ! Forward region of large grains: the diffraction peak and the first rings (period ~ pi/x)
      ! are much narrower than the angle grid. There S11, S12 are averaged over the angle bin of
      ! each grid angle, from >= 4 samples per ring period (at most 256 per grid interval); the
      ! number of bins treated this way covers theta < ~30/x. Elsewhere grid-point values.
      mfine = 1
      if (do_sca) then
        mfine = ceiling(4.0_r2 * sizepar * d_angle * deg2rad / pi)
        mfine = max(1, min(256, mfine))
        if (mfine > 1) mfine = 2 * ((mfine + 1) / 2)
      end if
      allocate(sf1(n_angle2), sf2(n_angle2), f11(n_angle2), f12(n_angle2))

      if (switch_grain_shape == 1) then
        ! distribution of hollow spheres: efficiencies and normalised S11, S12 averaged over f
        call dhs_observing(ri, sizepar, qext, qsca, qabs, f11, f12, ier)
        if (ier == 0) then
          if (.not. (qabs > 0.0_r2)) qabs = 0.0_r2
          Qabs_o(il, ig, is) = qabs
          Qsca_o(il, ig, is) = qsca
          Cabs(il, ig, is)   = qabs * pi * (rad(ig) * micron)**2
          Csca(il, ig, is)   = qsca * pi * (rad(ig) * micron)**2
          if (do_sca) then
            S11_norm(il, :, ig, is) = f11
            S12_norm(il, :, ig, is) = f12
          end if
          deallocate(sf1, sf2, f11, f12)
          cycle
        end if
      else
        call shexqnn2(ri, sizepar, qext, qsca, qabs, qbk, qpr, alb, gsc, ier, sf1, sf2, do_sca, n_angle, &
                      ab, nab)
      end if

      if (ier /= 0) then
        !$omp critical (mie_error)
        call dms_fail("Mie calculation failed at the observing wavelength "//real_to_string(wavelength_obs(il)))
        !$omp end critical (mie_error)
        deallocate(sf1, sf2, f11, f12)
        if (allocated(ab)) deallocate(ab)
        cycle
      end if
      if (.not. (qabs > 0.0_r2)) qabs = 0.0_r2

      Qabs_o(il, ig, is) = qabs
      Qsca_o(il, ig, is) = qsca
      Cabs(il, ig, is)   = qabs * pi * (rad(ig) * micron)**2
      Csca(il, ig, is)   = qsca * pi * (rad(ig) * micron)**2

      if (do_sca) then
        norm = pi * sizepar**2 * qsca                          ! k^2 Csca
        f11(:) = 0.5_r2 * (abs(sf2(:))**2 + abs(sf1(:))**2) / norm
        f12(:) = 0.5_r2 * (abs(sf2(:))**2 - abs(sf1(:))**2) / norm
        S11_norm(il, :, ig, is) = f11
        S12_norm(il, :, ig, is) = f12
        if (mfine > 1) then
          kfwd = min(n_angle2 - 1, max(2, ceiling(30.0_r2 / (sizepar * d_angle * deg2rad))))
          call forward_bins(ab, nab, mfine, kfwd, norm, S11_norm(il, :, ig, is), S12_norm(il, :, ig, is))
        end if
      end if
      deallocate(sf1, sf2, f11, f12)
      if (allocated(ab)) deallocate(ab)

    end do
    !$omp end parallel do

    deallocate(Qabs_o, Qsca_o)
    if (dms_ierr /= 0) return

    if (switch_miex_save == 1) call save_optics_cache()

  end if

  ! Each size bin is represented by one grain size. A coarse size grid leaves radial "spokes" in
  ! Mie Q_phi images and ripples with the period of the size grid in the (sub)millimetre SED
  if (n_grain_radius > 1 .and. grain_radius_max > grain_radius_min .and. verbose) then
    if (real(n_grain_radius, kind=r2) / log10(grain_radius_max / grain_radius_min) < 25.0_r2) then
      print *, "    note: ", real_to_string(real(n_grain_radius, kind=r2) / log10(grain_radius_max / &
               grain_radius_min)), " size bins per decade; use >= 30 for smooth Mie polarisation", &
               " images and (sub)millimetre fluxes"
    end if
  end if

  if (.not. do_sca) return

  !========================================================================================
  ! 2. Henyey-Greenstein phase functions
  !========================================================================================

  if (switch_phase_function >= 2) then

    do ia = 1, n_angle2
      ct = cos(x_angle(ia) * deg2rad)
      if (switch_phase_function == 2) then
        ph = hg_phase(hg_g1, ct)
      else
        ph = hg_w * hg_phase(hg_g1, ct) + (1.0_r2 - hg_w) * hg_phase(hg_g2, ct)
      end if
      pp = -pol_max * (1.0_r2 - ct*ct) / (1.0_r2 + ct*ct)
      S11_norm(:, ia, :, :) = ph
      S12_norm(:, ia, :, :) = pp * ph
    end do

    if (verbose) print *, "    Henyey-Greenstein phase function used for the scattered light"

  end if

  !========================================================================================
  ! 2b. energy-conserving phase function on the angle grid
  !     radiation_point interpolates S11 linearly in theta between the grid angles. Where the
  !     forward diffraction peak (width ~ 1/x rad) is not resolved by the grid, the linearly
  !     interpolated phase function does not contain exactly the scattered energy. S11 and S12
  !     are therefore rescaled so that int S11_lin dOmega = 1 (exact weights of the linear
  !     interpolation), but only in the forward region theta < 6/x (at least the first three grid
  !     angles), where the missing energy lies. For small grains and Henyey-Greenstein functions
  !     this region is the whole sphere: a uniform renormalisation that removes the quadrature
  !     error of the grid (~1e-5).
  !     (Earlier drafts put the whole difference into the value at theta = 0; for smooth
  !     phase functions this changed S11(0) by up to ~50 %.)
  !========================================================================================

  call phase_function_energy_conservation()

  !========================================================================================
  ! 3. scattering kernels
  !========================================================================================

  allocate(KS11(n_wavelength_obs, n_angle2, n_grain_radius, n_species), &
           KS12(n_wavelength_obs, n_angle2, n_grain_radius, n_species))

  do is = 1, n_species
    do ig = 1, n_grain_radius
      do ia = 1, n_angle2
        KS11(:, ia, ig, is) = Csca(:, ig, is) * S11_norm(:, ia, ig, is) * pi * Rstar**2 * star_I_obs(:)
        KS12(:, ia, ig, is) = Csca(:, ig, is) * S12_norm(:, ia, ig, is) * pi * Rstar**2 * star_I_obs(:)
      end do
    end do
  end do

  deallocate(S11_norm, S12_norm)

contains

  !---------------------------------------------------------------------------------------------
  ! DHS (Min et al. 2005) at one observing wavelength and grain size: for every vacuum fraction f
  ! the hollow sphere is treated like a Mie grain (amplitudes on the angle grid, bin averages in
  ! the forward region, normalisation); S11, S12 are then averaged with weights w_j C_sca,j, the
  ! efficiencies with w_j (relative to the volume-equivalent radius).
  !---------------------------------------------------------------------------------------------
  subroutine dhs_observing(ri, x_eq, qe, qs, qa, s11, s12, ier)

    complex(r2), intent(in)                 :: ri
    real(r2), intent(in)                    :: x_eq
    real(r2), intent(out)                   :: qe, qs, qa
    real(r2), dimension(:), intent(out)     :: s11, s12
    integer, intent(out)                    :: ier

    real(r2), dimension(dhs_nf) :: f, w
    complex(r2), allocatable, dimension(:,:) :: abj
    complex(r2), allocatable, dimension(:)   :: a1, a2
    real(r2), allocatable, dimension(:)      :: g11, g12
    real(r2) :: x_in, x_out, qej, qsj, qbj, qpj, gj, sw, nrm
    integer  :: j, nabj, mf, kf

    qe = 0.0_r2;  qs = 0.0_r2;  qa = 0.0_r2;  s11 = 0.0_r2;  s12 = 0.0_r2;  ier = 0
    allocate(a1(n_angle2), a2(n_angle2), g11(n_angle2), g12(n_angle2))
    call dhs_points(dhs_fmax, dhs_nf, f, w)
    do j = 1, dhs_nf
      call dhs_radii(x_eq, f(j), x_in, x_out)
      call coated_ab((1.0_r2, 0.0_r2), ri, x_in, x_out, abj, nabj, ier)
      if (ier /= 0) return
      call q_from_ab(abj, nabj, x_out, qej, qsj, qbj, qpj, gj)
      sw = w(j) * (x_out / x_eq)**2
      qe = qe + sw * qej
      qs = qs + sw * qsj
      if (do_sca .and. qsj > 0.0_r2) then
        call mie_amplitudes_grid(abj, nabj, n_angle, a1, a2)
        nrm = pi * x_out**2 * qsj                              ! k^2 Csca of this sphere
        g11 = 0.5_r2 * (abs(a2)**2 + abs(a1)**2) / nrm
        g12 = 0.5_r2 * (abs(a2)**2 - abs(a1)**2) / nrm
        mf = max(1, min(256, ceiling(4.0_r2 * x_out * d_angle * deg2rad / pi)))
        if (mf > 1) then
          mf = 2 * ((mf + 1) / 2)
          kf = min(n_angle2 - 1, max(2, ceiling(30.0_r2 / (x_out * d_angle * deg2rad))))
          call forward_bins(abj, nabj, mf, kf, nrm, g11, g12)
        end if
        s11 = s11 + sw * qsj * g11
        s12 = s12 + sw * qsj * g12
      end if
      deallocate(abj)
    end do
    qa = qe - qs
    if (qs > 0.0_r2) then
      s11 = s11 / qs
      s12 = s12 / qs
    end if
    deallocate(a1, a2, g11, g12)

  end subroutine dhs_observing

  !---------------------------------------------------------------------------------------------
  ! bin averages of S11, S12 for the grid angles j = 1 .. kf+1 (forward region):
  !   <S> = int S sin(theta) dtheta / int sin(theta) dtheta over [theta_j - d/2, theta_j + d/2]
  ! from mf samples per grid interval (trapezoidal rule), amplitudes from the Mie coefficients.
  !---------------------------------------------------------------------------------------------
  subroutine forward_bins(abc, nabc, mf, kf, nrm, s11, s12)

    complex(r2), dimension(:,:), intent(in) :: abc
    integer, intent(in)                     :: nabc, mf, kf
    real(r2), intent(in)                    :: nrm
    real(r2), dimension(:), intent(inout)   :: s11, s12
    integer  :: j, k, nf, k1, k2
    real(r2) :: dth, w, t, n11, n12, den
    real(r2),    allocatable :: muf(:), g11(:), g12(:)
    complex(r2), allocatable :: a1(:), a2(:)

    dth = d_angle * deg2rad / real(mf, kind=r2)
    nf  = kf * mf + mf / 2 + 1                       ! samples from 0 to (kf + 1/2) d
    allocate(muf(nf), g11(nf), g12(nf), a1(nf), a2(nf))
    do k = 1, nf
      muf(k) = cos(real(k - 1, kind=r2) * dth)
    end do
    call mie_amplitudes(abc, nabc, muf, a1, a2)
    g11 = 0.5_r2 * (abs(a2)**2 + abs(a1)**2) / nrm
    g12 = 0.5_r2 * (abs(a2)**2 - abs(a1)**2) / nrm

    do j = 1, kf + 1
      k1 = max(1,  (j - 1) * mf + 1 - mf / 2)
      k2 = min(nf, (j - 1) * mf + 1 + mf / 2)
      n11 = 0.0_r2;  n12 = 0.0_r2;  den = 0.0_r2
      do k = k1, k2
        w = dth
        if (k == k1 .or. k == k2) w = 0.5_r2 * dth
        t = real(k - 1, kind=r2) * dth
        n11 = n11 + w * g11(k) * sin(t)
        n12 = n12 + w * g12(k) * sin(t)
        den = den + w * sin(t)
      end do
      if (den > 0.0_r2) then
        s11(j) = n11 / den;  s12(j) = n12 / den
      end if
    end do
    deallocate(muf, g11, g12, a1, a2)

  end subroutine forward_bins

  subroutine phase_function_energy_conservation()

    real(r2), dimension(n_angle2) :: wp              ! int S_lin dOmega = sum wp(ia) S(ia)
    real(r2) :: ta, tb, dt, wa, wb, t_in, t_out, fac, x, dmax
    integer  :: jl, jg, js, ja, kf, nfwd

    wp = 0.0_r2
    do ja = 1, n_angle2 - 1
      ta = x_angle(ja) * deg2rad
      tb = x_angle(ja+1) * deg2rad
      dt = tb - ta
      wb = 2.0_r2 * pi * (sin(tb) - sin(ta) - dt * cos(tb)) / dt
      wa = 2.0_r2 * pi * (cos(ta) - cos(tb)) - wb
      wp(ja)   = wp(ja)   + wa
      wp(ja+1) = wp(ja+1) + wb
    end do

    dmax = 0.0_r2;  nfwd = 0
    do js = 1, n_species
      do jg = 1, n_grain_radius
        do jl = 1, n_wavelength_obs
          kf = n_angle2 - 1
          if (switch_phase_function == 1) then
            x  = 2.0_r2 * pi * rad(jg) / wavelength_obs(jl)
            kf = min(n_angle2 - 1, max(2, ceiling(6.0_r2 / (x * d_angle * deg2rad))))
          end if
          t_in  = sum(wp(1:kf+1) * S11_norm(jl, 1:kf+1, jg, js))
          t_out = sum(wp(kf+2:n_angle2) * S11_norm(jl, kf+2:n_angle2, jg, js))
          if (t_in > 0.0_r2 .and. 1.0_r2 - t_out > 0.0_r2) then
            fac = (1.0_r2 - t_out) / t_in
            S11_norm(jl, 1:kf+1, jg, js) = S11_norm(jl, 1:kf+1, jg, js) * fac
            S12_norm(jl, 1:kf+1, jg, js) = S12_norm(jl, 1:kf+1, jg, js) * fac
          else if (t_in + t_out > 0.0_r2) then
            fac = 1.0_r2 / (t_in + t_out)
            S11_norm(jl, :, jg, js) = S11_norm(jl, :, jg, js) * fac
            S12_norm(jl, :, jg, js) = S12_norm(jl, :, jg, js) * fac
          else
            fac = 1.0_r2
          end if
          if (kf < n_angle2 - 1) nfwd = nfwd + 1
          dmax = max(dmax, abs(fac - 1.0_r2))
        end do
      end do
    end do

    if (verbose) then
      print *, "    phase function normalised on the angle grid (largest correction ", &
               real_to_string(100.0_r2 * dmax), " %; forward region only for ", nfwd, " of ", &
               n_wavelength_obs * n_grain_radius * n_species, " (lambda, a, species))"
    end if

  end subroutine phase_function_energy_conservation

end subroutine Qabs_Qsca_S11_S12_cals
