subroutine particle_abundances

  !---------------------------------------------------------------------------------------------------
  ! Number of grains per unit (initial) dust mass at every distance, including
  !   a) sublimation  (grains whose equilibrium temperature >= T_subl of their component either
  !                    vanish or become grains of another component; the components are processed
  !                    in order of increasing T_subl so that a residue can itself sublimate later)
  !   b) blow-out     (grains with beta = F_rad / F_grav >= beta_crit are removed; beta_crit = 1/2
  !                    for release from circular orbits, (1-e)/2 or (1+e)/2 for release at the
  !                    pericentre or apocentre of an eccentric belt, or a value given in the input),
  ! following the scheme of the debris-disk version with sublimation (Kirchschlager & Wolf 2013 for
  ! beta). The initial size distribution, normalised to 1 kg of dust:
  !     npm0(ig, ic) = abun(ic) w(ig) / sum_{ig,ic} abun w m_grain
  !
  ! Kernels per kg of initial dust mass at every distance id (used by radiation_point):
  !   j_reem(il, id)  = sum N Cabs(lambda) B_lambda(T)                       [W/m/sr/kg]
  !   fl_dist(id)     = sum N pi a^2 <Qabs>*      (fractional luminosity)    [m^2/kg]
  !   sig_dist(id)    = sum N pi a^2              (geometrical cross section)[m^2/kg]
  !   mfrac_dist(id)  = surviving mass fraction
  !   Ksca_sum, Kpol_sum: distance-independent scattering kernels if the grain population does not
  !                       depend on distance (no sublimation, no T-dependent optics)
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2), allocatable, dimension(:,:)   :: npm0, beta_sp
  real(r2), allocatable, dimension(:)     :: Istar_ri
  integer,  allocatable, dimension(:)     :: order
  integer  :: ic, ig, id, is, k, kk, itmp, ic_from, ic_to, il, ia, nremoved
  real(r2) :: m0, Tg, w
  character(len=:), allocatable :: crit_txt

  if (verbose) print *, " "
  if (verbose) print *, " 6) Grain populations (sublimation, blow-out) and emission kernels"
  if (verbose) print *, " "

  !------------------------------------------------------------------------------------------------
  ! initial population per kg
  !------------------------------------------------------------------------------------------------

  allocate(npm0(n_grain_radius, n_grain_components))
  m0 = 0.0_r2
  do ic = 1, n_grain_components
    npm0(:, ic) = abun(ic) * size_weight(:)
    m0 = m0 + sum(npm0(:, ic) * grain_mass(:, ic))
  end do
  npm0(:, :) = npm0(:, :) / m0

  allocate(npm(n_grain_radius, n_grain_components, n_dist))
  do id = 1, n_dist
    npm(:, :, id) = npm0(:, :)
  end do

  !------------------------------------------------------------------------------------------------
  ! a) sublimation
  !------------------------------------------------------------------------------------------------

  if (switch_sublimation == 1) then

    do ic = 1, n_grain_components
      if (subl_to_comp(ic) /= 0 .and. T_subl(ic) > 0.0_r2) then
        if (subl_to_comp(ic) == ic) then
          call dms_fail("component "//integer_to_string(ic)//" cannot sublimate into itself")
          return
        end if
        if (T_subl(subl_to_comp(ic)) > 0.0_r2 .and. T_subl(subl_to_comp(ic)) <= T_subl(ic)) then
          call dms_fail("component "//integer_to_string(ic)//" sublimates into component "// &
                        integer_to_string(subl_to_comp(ic))//" whose sublimation temperature is not higher")
          return
        end if
      end if
    end do

    ! order of increasing sublimation temperature
    allocate(order(n_grain_components))
    order = [(ic, ic = 1, n_grain_components)]
    do k = 1, n_grain_components - 1
      do kk = k + 1, n_grain_components
        if (sort_key(order(kk)) < sort_key(order(k))) then
          itmp = order(k); order(k) = order(kk); order(kk) = itmp
        end if
      end do
    end do

    do k = 1, n_grain_components
      ic_from = order(k)
      if (T_subl(ic_from) <= 0.0_r2) cycle
      ic_to = subl_to_comp(ic_from)
      do id = 1, n_dist
        do ig = 1, n_grain_radius
          if (npm(ig, ic_from, id) <= 0.0_r2) cycle
          if (T_dist(ig, ic_from, id) >= T_subl(ic_from)) then
            if (ic_to > 0) npm(ig, ic_to, id) = npm(ig, ic_to, id) + npm(ig, ic_from, id)
            npm(ig, ic_from, id) = 0.0_r2
          end if
        end do
      end do
    end do

    do ic = 1, n_grain_components
      if (T_subl(ic) <= 0.0_r2) cycle
      kk = 0
      do id = n_dist, 1, -1
        if (any(T_dist(:, ic, id) >= T_subl(ic))) then
          kk = id
          exit
        end if
      end do
      if (kk == 0) then
        if (verbose) print *, "    component ", ic, ": no sublimation (T < T_subl everywhere)"
      else if (kk == n_dist) then
        if (verbose) print *, "    component ", ic, ": (some) grains sublimate everywhere in the disk"
      else
        if (verbose) print *, "    component ", ic, ": sublimation zone inside ", &
                 real_to_string(dist_grid(kk) / au), " au (T_subl = ", real_to_string(T_subl(ic)), " K)"
      end if
    end do

  end if

  !------------------------------------------------------------------------------------------------
  ! b) radiation pressure: beta = 3 pi R*^2 int Qpr I* dlambda / (4 G M* c rho a)
  !------------------------------------------------------------------------------------------------

  ! beta for every optical species (temperature-dependent optical constants: the species used at a
  ! given distance decides); reported beta_rp(ig, ic) = largest value over the distances
  allocate(beta_rp(n_grain_radius, n_grain_components), beta_sp(n_grain_radius, n_species))
  allocate(Istar_ri(n_wavelength_int))
  do k = 1, n_wavelength_int
    Istar_ri(k) = stellar_radiance(wavelength_int(k) * micron)
  end do

  do is = 1, n_species
    ic = species_comp(is)
    do ig = 1, n_grain_radius
      if (Mstar > 0.0_r2) then
        beta_sp(ig, is) = 3.0_r2 * pi * Rstar**2 * sum(Qpr_ri(:, ig, is) * Istar_ri(:) * dlambda_int(:)) / &
                          (4.0_r2 * G_grav * c * Mstar * particle_bulk_density(ic) * 1000.0_r2 &
                           * rad(ig) * micron)
      else
        beta_sp(ig, is) = 0.0_r2
      end if
    end do
  end do

  do ic = 1, n_grain_components
    do ig = 1, n_grain_radius
      beta_rp(ig, ic) = 0.0_r2
      do id = 1, n_dist
        beta_rp(ig, ic) = max(beta_rp(ig, ic), beta_sp(ig, sp_dist(ig, ic, id)))
      end do
    end do
  end do

  call blowout_thresholds

  if (switch_blowout == 1) then
    if (Mstar <= 0.0_r2) then
      call dms_fail("blow-out requires a positive stellar mass")
      return
    end if
    if (beta_crit_max - beta_crit_min > 1.0e-12_r2) then
      crit_txt = real_to_string(beta_crit_min)//" - "//real_to_string(beta_crit_max)
    else
      crit_txt = real_to_string(beta_crit_min)
    end if
    if (verbose .and. (switch_blowout_release /= 1 .or. blowout_beta_user > 0.0_r2)) &
      print *, "    blow-out threshold beta_crit = ", crit_txt
    do ic = 1, n_grain_components
      nremoved = 0
      do ig = 1, n_grain_radius
        do id = 1, n_dist
          if (beta_sp(ig, sp_dist(ig, ic, id)) >= beta_crit_dist(id)) npm(ig, ic, id) = 0.0_r2
        end do
        if (beta_rp(ig, ic) >= beta_crit_min) nremoved = nremoved + 1
      end do
      if (verbose) print *, "    component ", ic, ": ", nremoved, " size bins blown out (beta >= ", crit_txt, ")"
    end do
  else
    do ic = 1, n_grain_components
      if (any(beta_rp(:, ic) >= 0.5_r2 .and. npm0(:, ic) > 0.0_r2)) then
        if (verbose) print *, "    note: component ", ic, " contains grains with beta >= 0.5 (not removed)"
      end if
    end do
  end if

  if (all(npm <= 0.0_r2)) then
    if (switch_blowout == 1 .and. all(beta_rp >= beta_crit_max .or. npm0 <= 0.0_r2)) then
      call dms_fail("no grains left: all grain sizes are blown out (beta >= "//real_to_string(beta_crit_max)// &
                    "); increase the grain sizes")
    else
      call dms_fail("no grains left after sublimation / blow-out at any distance")
    end if
    return
  end if

  ! one population for all distances only without sublimation, T-dependent optics and a
  ! distance-dependent blow-out threshold
  uniform_population = (switch_sublimation == 0 .and. n_species == n_grain_components)
  if (switch_blowout == 1 .and. beta_crit_max - beta_crit_min > 1.0e-12_r2) uniform_population = .false.

  !------------------------------------------------------------------------------------------------
  ! kernels per kg of initial dust
  !------------------------------------------------------------------------------------------------

  allocate(j_reem(n_wavelength_obs, n_dist), fl_dist(n_dist), sig_dist(n_dist), mfrac_dist(n_dist))

  !$omp parallel do schedule(static) private(id, ic, ig, is, Tg, w, il)
  do id = 1, n_dist
    j_reem(:, id)  = 0.0_r2
    fl_dist(id)    = 0.0_r2
    sig_dist(id)   = 0.0_r2
    mfrac_dist(id) = 0.0_r2
    do ic = 1, n_grain_components
      do ig = 1, n_grain_radius
        w = npm(ig, ic, id)
        if (w <= 0.0_r2) cycle
        is = sp_dist(ig, ic, id)
        Tg = T_dist(ig, ic, id)
        mfrac_dist(id) = mfrac_dist(id) + w * grain_mass(ig, ic)
        sig_dist(id)   = sig_dist(id)   + w * pi * (rad(ig) * micron)**2
        fl_dist(id)    = fl_dist(id)    + w * pi * (rad(ig) * micron)**2 * Qstar_mean(ig, is)
        if (do_reem) then
          do il = 1, n_wavelength_obs
            j_reem(il, id) = j_reem(il, id) + w * Cabs(il, ig, is) * Planck_fn(Tg, lambda_obs_m(il))
          end do
        end if
      end do
    end do
  end do
  !$omp end parallel do

  if (do_sca .and. uniform_population) then
    allocate(Ksca_sum(n_wavelength_obs, n_angle2), Kpol_sum(n_wavelength_obs, n_angle2))
    Ksca_sum = 0.0_r2
    Kpol_sum = 0.0_r2
    do ic = 1, n_grain_components
      is = comp_first_species(ic)
      do ig = 1, n_grain_radius
        w = npm(ig, ic, 1)
        if (w <= 0.0_r2) cycle
        do ia = 1, n_angle2
          Ksca_sum(:, ia) = Ksca_sum(:, ia) + w * KS11(:, ia, ig, is)
          Kpol_sum(:, ia) = Kpol_sum(:, ia) + w * KS12(:, ia, ig, is)
        end do
      end do
    end do
  end if

  if (uniform_population) then
    if (verbose) print *, "    grain population independent of distance (fast scattering kernel)"
  else
    if (verbose) print *, "    grain population depends on distance (sublimation / T-dependent optics)"
  end if

  deallocate(npm0, beta_sp, Istar_ri)

contains

  !-------------------------------------------------------------------------------------------------
  ! Blow-out threshold at every distance of the distance grid. For release at the pericentre or
  ! apocentre the threshold of belt k is (1 -+ e_k)/2; where belts overlap in distance, the smallest
  ! threshold applies, and outside all belts (no dust) the smallest threshold of all belts.
  !-------------------------------------------------------------------------------------------------
  subroutine blowout_thresholds

    integer  :: kb
    real(r2) :: crit_belt(max_belts), crit_all, r_hi, s, cmin

    if (allocated(beta_crit_dist)) deallocate(beta_crit_dist)
    allocate(beta_crit_dist(n_dist))

    do kb = 1, n_belts
      if (blowout_beta_user > 0.0_r2) then
        crit_belt(kb) = blowout_beta_user
      else
        select case (switch_blowout_release)
        case (2);      crit_belt(kb) = 0.5_r2 * (1.0_r2 - belt_ecc(kb))
        case (3);      crit_belt(kb) = 0.5_r2 * (1.0_r2 + belt_ecc(kb))
        case default;  crit_belt(kb) = 0.5_r2
        end select
      end if
    end do
    crit_all = minval(crit_belt(1:n_belts))

    do id = 1, n_dist
      s    = dist_grid(id)
      cmin = huge(1.0_r2)
      do kb = 1, n_belts
        r_hi = belt_rout_eff(kb)
        if (switch_vertical == 2) r_hi = r_hi / cos(open_angle)
        if (s >= belt_rin_eff(kb) * (1.0_r2 - 1.0e-12_r2) .and. s <= r_hi * (1.0_r2 + 1.0e-12_r2)) &
          cmin = min(cmin, crit_belt(kb))
      end do
      if (cmin > 1.0e30_r2) cmin = crit_all
      beta_crit_dist(id) = cmin
    end do

    beta_crit_min = minval(beta_crit_dist)
    beta_crit_max = maxval(beta_crit_dist)

  end subroutine blowout_thresholds


  pure function sort_key(icomp) result(key)
    integer, intent(in) :: icomp
    real(r2)            :: key
    key = T_subl(icomp)
    if (key <= 0.0_r2) key = huge(1.0_r2)
  end function sort_key

end subroutine particle_abundances
