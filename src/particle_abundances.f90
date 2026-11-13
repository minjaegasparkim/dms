subroutine particle_abundances

  !---------------------------------------------------------------------------------------------------
  ! Number of grains per unit (initial) dust mass at every distance, including
  !   0) photosputtering (dms 2.2; grains of an icy component that lose their ice before they are
  !                    destroyed by collisions, drift inwards or reach the age of the system become
  !                    grains of another component or vanish; Grigorieva et al. 2007)
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
  real(r2) :: m0, Tg, w, Tlo, Thi
  character(len=:), allocatable :: crit_txt

  if (verbose) print *, " "
  if (verbose) print *, " 6) Grain populations (photosputtering, sublimation, blow-out) and emission kernels"
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
  ! radiation pressure: beta = 3 pi R*^2 int Qpr I* dlambda / (4 G M* c rho a)
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

  !------------------------------------------------------------------------------------------------
  ! photosputtering (dms 2.2): icy grains that lose their ice before the reference time become
  ! grains of the residue component; done before sublimation so that a residue can sublimate
  !------------------------------------------------------------------------------------------------

  orbit_mode = .false.;  n_obins = 0;  n_ocell = 0;  ob_first = 0;  ob_count = 0
  if (switch_photosputtering == 1) then
    call blowout_thresholds
    if (switch_ps_orbit == 2 .and. any(belt_ecc(1:n_belts) > 0.0_r2)) then
      call setup_orbit_bins
      if (dms_ierr /= 0) return
    end if
    call photosputtering
    if (dms_ierr /= 0) return
  end if

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
  ! b) blow-out
  !------------------------------------------------------------------------------------------------

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

  !------------------------------------------------------------------------------------------------
  ! grains colder or hotter than the data sets of their component (temperature-dependent data):
  ! the coldest or hottest set is used; reported once per component
  !------------------------------------------------------------------------------------------------

  if (allocated(comp_T_outside)) deallocate(comp_T_outside)
  allocate(comp_T_outside(4, n_grain_components))
  comp_T_outside = -1.0_r2
  do ic = 1, n_grain_components
    if (comp_n_sets(ic) < 2) cycle
    Tlo = set_T(comp_first_set(ic))
    Thi = set_T(comp_first_set(ic) + comp_n_sets(ic) - 1)
    do id = 1, n_dist
      do ig = 1, n_grain_radius
        if (npm(ig, ic, id) <= 0.0_r2) cycle
        Tg = T_dist(ig, ic, id)
        if (Tg < Tlo) then
          if (comp_T_outside(1, ic) < 0.0_r2) then
            comp_T_outside(1:2, ic) = Tg
          else
            comp_T_outside(1, ic) = min(comp_T_outside(1, ic), Tg)
            comp_T_outside(2, ic) = max(comp_T_outside(2, ic), Tg)
          end if
        else if (Tg > Thi) then
          if (comp_T_outside(3, ic) < 0.0_r2) then
            comp_T_outside(3:4, ic) = Tg
          else
            comp_T_outside(3, ic) = min(comp_T_outside(3, ic), Tg)
            comp_T_outside(4, ic) = max(comp_T_outside(4, ic), Tg)
          end if
        end if
      end do
    end do
    if (verbose .and. comp_T_outside(1, ic) > 0.0_r2) &
      print *, "    warning: component ", ic, ": grain temperatures ", real_to_string(comp_T_outside(1, ic)), &
               " - ", real_to_string(comp_T_outside(2, ic)), " K are below the coldest data set (", &
               real_to_string(Tlo), " K); that set is used"
    if (verbose .and. comp_T_outside(3, ic) > 0.0_r2) &
      print *, "    warning: component ", ic, ": grain temperatures ", real_to_string(comp_T_outside(3, ic)), &
               " - ", real_to_string(comp_T_outside(4, ic)), " K are above the hottest data set (", &
               real_to_string(Thi), " K); that set is used"
  end do

  ! one population for all distances only without sublimation, photosputtering, T-dependent optics and a
  ! distance-dependent blow-out threshold
  uniform_population = (switch_sublimation == 0 .and. switch_photosputtering == 0 .and. &
                        n_species == n_grain_components)
  if (switch_blowout == 1 .and. beta_crit_max - beta_crit_min > 1.0e-12_r2) uniform_population = .false.

  !------------------------------------------------------------------------------------------------
  ! kernels per kg of initial dust
  !------------------------------------------------------------------------------------------------

  ! distance entries 1..n_dist, then the cells of the orbit bins (orbit_populations.f90)
  allocate(j_reem(n_wavelength_obs, n_dist + n_ocell), fl_dist(n_dist + n_ocell), &
           sig_dist(n_dist + n_ocell), mfrac_dist(n_dist + n_ocell))

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

  if (orbit_mode) then
    call build_orbit_cells
    if (dms_ierr /= 0) return
  end if

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
    if (verbose) print *, "    grain population depends on distance (sublimation / photosputtering / T-dependent optics)"
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
    real(r2) :: crit_belt(max_belts), crit_all, r_hi, s, cmin, e_all, e_here

    if (allocated(beta_crit_dist)) deallocate(beta_crit_dist)
    if (allocated(ecc_parent_dist)) deallocate(ecc_parent_dist)
    allocate(beta_crit_dist(n_dist), ecc_parent_dist(n_dist))

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
    e_all = belt_ecc(minloc(crit_belt(1:n_belts), dim=1))

    do id = 1, n_dist
      s    = dist_grid(id)
      cmin = huge(1.0_r2)
      e_here = e_all
      do kb = 1, n_belts
        r_hi = belt_rout_eff(kb)
        if (switch_vertical == 2) r_hi = r_hi / cos(open_angle)
        if (s >= belt_rin_eff(kb) * (1.0_r2 - 1.0e-12_r2) .and. s <= r_hi * (1.0_r2 + 1.0e-12_r2)) then
          if (crit_belt(kb) < cmin) e_here = belt_ecc(kb)
          cmin = min(cmin, crit_belt(kb))
        end if
      end do
      if (cmin > 1.0e30_r2) cmin = crit_all
      beta_crit_dist(id) = cmin
      ! parent orbits there (used by photosputtering with ps_orbit = 'release'): the belt that sets
      ! the threshold; circular release: e = 0
      ecc_parent_dist(id) = merge(0.0_r2, e_here, switch_blowout_release == 1)
    end do

    beta_crit_min = minval(beta_crit_dist)
    beta_crit_max = maxval(beta_crit_dist)

  end subroutine blowout_thresholds


  !-------------------------------------------------------------------------------------------------
  ! Orbit bins of the eccentric belts (photosputtering with ps_orbit = 'release'): ps_orbit_bins
  ! logarithmic bins in semi-major axis between belt_Rmin and belt_Rmax; every bin keeps the
  ! distance indices its orbits can reach (a (1 - e) to a (1 + e), within the belt edges)
  !-------------------------------------------------------------------------------------------------
  subroutine setup_orbit_bins

    integer  :: kb, j, ob, nb_k, off
    real(r2) :: s_lo, s_hi, gb

    n_obins = 0
    do kb = 1, n_belts
      if (belt_ecc(kb) <= 0.0_r2) cycle
      nb_k = max(1, ps_nbins)
      if (belt_Rmax(kb) <= belt_Rmin(kb) * (1.0_r2 + 1.0e-12_r2)) nb_k = 1
      ob_first(kb) = n_obins + 1
      ob_count(kb) = nb_k
      n_obins = n_obins + nb_k
    end do
    if (n_obins == 0) return

    if (allocated(ob_belt)) deallocate(ob_belt, ob_idlo, ob_idhi, ob_off, ob_alo, ob_ahi)
    allocate(ob_belt(n_obins), ob_idlo(n_obins), ob_idhi(n_obins), ob_off(n_obins), &
             ob_alo(n_obins), ob_ahi(n_obins))
    off = 0
    do kb = 1, n_belts
      if (ob_count(kb) <= 0) cycle
      do j = 1, ob_count(kb)
        ob = ob_first(kb) + j - 1
        ob_belt(ob) = kb
        ob_alo(ob) = belt_Rmin(kb) * (belt_Rmax(kb) / belt_Rmin(kb))**(real(j - 1, kind=r2) / real(ob_count(kb), kind=r2))
        ob_ahi(ob) = belt_Rmin(kb) * (belt_Rmax(kb) / belt_Rmin(kb))**(real(j, kind=r2) / real(ob_count(kb), kind=r2))
        s_lo = max(ob_alo(ob) * (1.0_r2 - belt_ecc(kb)), belt_rin_eff(kb))
        s_hi = min(ob_ahi(ob) * (1.0_r2 + belt_ecc(kb)), belt_rout_eff(kb))
        if (switch_vertical == 2) then
          s_hi = s_hi / cos(open_angle)
        end if
        ob_idlo(ob) = max(1, dist_index(s_lo) - 1)
        ob_idhi(ob) = min(n_dist, max(ob_idlo(ob), dist_index(s_hi) + 1))
        ob_off(ob)  = off
        off = off + ob_idhi(ob) - ob_idlo(ob) + 1
      end do
    end do
    n_ocell = off

    gb = real(n_ocell, kind=r2) * real(n_wavelength_obs + 4, kind=r2) * 8.0_r2 / 1.0e9_r2
    if (gb > 8.0_r2) then
      call dms_fail("orbit bins of the eccentric belts need "//real_to_string(gb)// &
                    " GB; reduce ps_orbit_bins or the number of wavelengths")
      return
    end if
    if (allocated(ps_lost_ob)) deallocate(ps_lost_ob)
    allocate(ps_lost_ob(n_grain_radius, n_grain_components, n_obins))
    ps_lost_ob = .false.
    orbit_mode = .true.
    if (verbose) print *, "    photosputtering on release orbits: ", n_obins, " orbit bins in the eccentric ", &
                          "belt(s), ", n_ocell, " cells (", real_to_string(gb), " GB)"

  end subroutine setup_orbit_bins


  !-------------------------------------------------------------------------------------------------
  ! Final component of every initial grain at every distance, with and without ice loss
  ! (sublimation chain, then blow-out, as for the distance entries), and the kernels of the cells
  !-------------------------------------------------------------------------------------------------
  subroutine build_orbit_cells

    integer  :: jd, jg, src, cc, ob, kx, js, jl, jc
    real(r2) :: w, Tj

    if (allocated(dest_keep)) deallocate(dest_keep, dest_lost)
    allocate(dest_keep(n_grain_radius, n_grain_components, n_dist), &
             dest_lost(n_grain_radius, n_grain_components, n_dist))
    if (allocated(npm_init)) deallocate(npm_init)
    allocate(npm_init(n_grain_radius, n_grain_components))
    npm_init(:, :) = npm0(:, :)

    !$omp parallel do schedule(static) private(jd, jg, src)
    do jd = 1, n_dist
      do src = 1, n_grain_components
        do jg = 1, n_grain_radius
          dest_keep(jg, src, jd) = final_component(jg, src, jd)
          if (ps_Y(src) > 0.0_r2) then
            dest_lost(jg, src, jd) = 0
            if (ps_to_comp(src) > 0) dest_lost(jg, src, jd) = final_component(jg, ps_to_comp(src), jd)
          else
            dest_lost(jg, src, jd) = dest_keep(jg, src, jd)
          end if
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do schedule(dynamic) private(ob, jd, kx, src, jg, cc, w, js, Tj, jl)
    do ob = 1, n_obins
      do jd = ob_idlo(ob), ob_idhi(ob)
        kx = n_dist + ob_off(ob) + jd - ob_idlo(ob) + 1
        j_reem(:, kx) = 0.0_r2
        fl_dist(kx) = 0.0_r2;  sig_dist(kx) = 0.0_r2;  mfrac_dist(kx) = 0.0_r2
        do src = 1, n_grain_components
          do jg = 1, n_grain_radius
            w = npm_init(jg, src)
            if (w <= 0.0_r2) cycle
            if (ps_lost_ob(jg, src, ob)) then
              cc = dest_lost(jg, src, jd)
            else
              cc = dest_keep(jg, src, jd)
            end if
            if (cc <= 0) cycle
            js = sp_dist(jg, cc, jd)
            Tj = T_dist(jg, cc, jd)
            mfrac_dist(kx) = mfrac_dist(kx) + w * grain_mass(jg, cc)
            sig_dist(kx)   = sig_dist(kx)   + w * pi * (rad(jg) * micron)**2
            fl_dist(kx)    = fl_dist(kx)    + w * pi * (rad(jg) * micron)**2 * Qstar_mean(jg, js)
            if (do_reem) then
              do jl = 1, n_wavelength_obs
                j_reem(jl, kx) = j_reem(jl, kx) + w * Cabs(jl, jg, js) * Planck_fn(Tj, lambda_obs_m(jl))
              end do
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    jc = 0
    do src = 1, n_grain_components
      if (ps_Y(src) > 0.0_r2) jc = jc + count(ps_lost_ob(:, src, :))
    end do
    if (verbose) print *, "    orbit bins: ", jc, " (size, orbit bin) pairs lose their ice"

  end subroutine build_orbit_cells


  ! component a grain of component c0 (size jg) ends as at distance jd: the sublimation chain
  ! (T >= T_subl -> sublimation_residue), then blow-out (0: removed)
  pure function final_component(jg, c0, jd) result(cc)
    integer, intent(in) :: jg, c0, jd
    integer             :: cc
    cc = c0
    if (switch_sublimation == 1) then
      do while (cc > 0)
        if (T_subl(cc) <= 0.0_r2) exit
        if (T_dist(jg, cc, jd) < T_subl(cc)) exit
        cc = subl_to_comp(cc)
      end do
    end if
    if (cc > 0 .and. switch_blowout == 1) then
      if (beta_sp(jg, sp_dist(jg, cc, jd)) >= beta_crit_dist(jd)) cc = 0
    end if
  end function final_component


  !-------------------------------------------------------------------------------------------------
  ! Ice loss by photosputtering (Grigorieva et al. 2007, A&A 475, 755). The ice of a grain of radius
  ! a is removed in
  !     t_ps = 4 rho_w a / (m_H2O Y N_abs),  N_abs = int_band Q(lambda) F_lambda(r) lambda/(h c) dlambda
  ! (the ice front recedes at m_H2O Y N_abs / (4 rho_w), their eq. 7 with eta = 1; rho_w = ice mass
  ! per grain volume; N_abs absorbed photons per unit cross section and time, plus the interstellar
  ! photons ps_isrf_flux times the band mean of Q). The grain loses its ice where t_ps < t_ref:
  !     collision: t_coll = [(s/s1)^-2 + (s/s2)^2.7] / (tau Omega) for s < s2, (s/s2)^0.3 / (tau Omega)
  !                above (their eqs 19-20), s1 = 1.2 s_1/2, s2 = 100 s_1/2, s_1/2 = largest radius
  !                with beta >= 1/2 at that distance (the smallest grain if there is none)
  !     pr:        t_PR = r^2 c / (2 beta G M*) (their eq. 16)
  !     age:       ps_age_yr
  !     shortest:  the shortest of these (collision only with ps_tau_v > 0, age only with ps_age_yr > 0)
  ! F_lambda: the stellar model, or the spectrum ps_uv_file where it covers the band.
  ! ps_orbit = 'release': the stellar photons are averaged over the Kepler orbit of the grain after
  ! its release at r from a parent orbit (a_p, e_p) with a_p = r, r/(1-e), r/(1+e) for release from
  ! a circular orbit, at pericentre, at apocentre (blowout_release, e of the belt):
  !     1/a_g = (1/a_p - 2 beta/r) / (1 - beta),   1 - e_g^2 = a_p (1 - e_p^2) / ((1 - beta) a_g),
  !     <r^-2> = 1 / (a_g^2 sqrt(1 - e_g^2))       (N_abs scaled by r^2 <r^-2>),
  ! and the P-R time is that of this orbit, a_g^2 c (1-e_g^2)^(3/2) / (beta G M* (2 + 3 e_g^2))
  ! (orbit-averaged decay of a, Wyatt & Whipple 1950). Unbound grains (1/a_g <= 0) are skipped. The
  ! collisional lifetime and the age are not changed; the grains stay at their release distance.
  ! With blow-out on, grains above the blow-out threshold are not tested (they leave as they are).
  !-------------------------------------------------------------------------------------------------
  subroutine photosputtering

    integer,  parameter :: nb = 400
    real(r2), parameter :: m_h2o = 18.01528_r2 * 1.66053906660e-27_r2        ! [kg]
    real(r2), parameter :: year  = 365.25_r2 * 86400.0_r2                    ! [s]
    real(r2) :: lam_b(nb), w_b(nb), phot_b(nb), q_b(nb), t_b(nb)
    integer  :: i_b(nb)
    real(r2), allocatable, dimension(:,:) :: tab, Gq, Qm
    real(r2), allocatable, dimension(:)   :: uv_lam, uv_F, loglam
    integer  :: kb, n_uv, ics, icto, idd, igg, iss, nshow, ks, ob
    real(r2) :: r, rho_w, s_h, x, fl, lam, a_p, e_p, a_c
    logical  :: use_coll, use_pr, use_age

    if (Mstar <= 0.0_r2) then
      call dms_fail("photosputtering requires a positive stellar mass");  return
    end if

    ! band grid (logarithmic) and trapezoidal weights [m]
    do kb = 1, nb
      lam_b(kb) = ps_lam1 * (ps_lam2 / ps_lam1)**(real(kb - 1, kind=r2) / real(nb - 1, kind=r2))
    end do
    do kb = 1, nb
      if (kb == 1) then
        w_b(kb) = 0.5_r2 * (lam_b(2) - lam_b(1)) * micron
      else if (kb == nb) then
        w_b(kb) = 0.5_r2 * (lam_b(nb) - lam_b(nb-1)) * micron
      else
        w_b(kb) = 0.5_r2 * (lam_b(kb+1) - lam_b(kb-1)) * micron
      end if
    end do

    ! UV spectrum (micron, F_nu [Jy] at the distance of the star), converted to F_lambda at 1 au
    n_uv = 0
    if (len_trim(ps_uv_name) > 0) then
      call read_table(trim(ps_uv_name), 2, tab, n_uv, 2)
      if (dms_ierr /= 0) return
      allocate(uv_lam(n_uv), uv_F(n_uv))
      uv_lam(:) = tab(1, :)
      if (any(uv_lam <= 0.0_r2) .or. any(uv_lam(2:n_uv) <= uv_lam(1:n_uv-1))) then
        call dms_fail("ps_uv_file: the wavelengths must be positive and increasing");  return
      end if
      if (any(tab(2, :) < 0.0_r2)) then
        call dms_fail("ps_uv_file: negative flux");  return
      end if
      uv_F(:) = tab(2, :) * 1.0e-26_r2 * c / (uv_lam(:) * micron)**2 * (dist_star / au)**2
      deallocate(tab)
    end if

    ! photon spectral flux at 1 au [photons / m^2 / s / m]
    ps_n_model = 0
    do kb = 1, nb
      lam = lam_b(kb) * micron
      if (n_uv > 0 .and. lam_b(kb) >= uv_lam(1) * (1.0_r2 - 1.0e-9_r2) .and. &
          lam_b(kb) <= uv_lam(n_uv) * (1.0_r2 + 1.0e-9_r2)) then
        fl = interp_lin(log(uv_lam), uv_F, log(lam_b(kb)))
      else
        if (n_uv > 0) ps_n_model = ps_n_model + 1
        fl = pi * stellar_radiance(lam) * (Rstar / au)**2
      end if
      phot_b(kb) = fl * lam / (h * c)
    end do
    ps_photons_1au = sum(w_b(:) * phot_b(:))

    ! absorption efficiency on the band grid (interpolation weights into the optical-data grid)
    ps_q_extrap = (ps_lam1 < wavelength_int(1) * (1.0_r2 - 1.0e-9_r2) .or. &
                   ps_lam2 > wavelength_int(n_wavelength_int) * (1.0_r2 + 1.0e-9_r2))
    if (switch_ps_absorption == 1 .and. n_wavelength_int >= 2) then
      allocate(loglam(n_wavelength_int))
      loglam(:) = log(wavelength_int(:))
      do kb = 1, nb
        x = log(min(max(lam_b(kb), wavelength_int(1)), wavelength_int(n_wavelength_int)))
        i_b(kb) = min(locate(loglam, x), n_wavelength_int - 1)
        t_b(kb) = (x - loglam(i_b(kb))) / (loglam(i_b(kb) + 1) - loglam(i_b(kb)))
        t_b(kb) = max(0.0_r2, min(1.0_r2, t_b(kb)))
      end do
    end if
    allocate(Gq(n_grain_radius, n_species), Qm(n_grain_radius, n_species))
    do iss = 1, n_species
      do igg = 1, n_grain_radius
        if (switch_ps_absorption == 2) then
          q_b(:) = 1.0_r2
        else if (n_wavelength_int < 2) then
          q_b(:) = max(Qabs_ri(1, igg, iss), 0.0_r2)
        else
          do kb = 1, nb
            q_b(kb) = max((1.0_r2 - t_b(kb)) * Qabs_ri(i_b(kb), igg, iss) + &
                          t_b(kb) * Qabs_ri(i_b(kb) + 1, igg, iss), 0.0_r2)
          end do
        end if
        Gq(igg, iss) = sum(w_b(:) * phot_b(:) * q_b(:))
        if (ps_photons_1au > 0.0_r2) then
          Qm(igg, iss) = Gq(igg, iss) / ps_photons_1au
        else
          Qm(igg, iss) = sum(q_b) / real(nb, kind=r2)
        end if
      end do
    end do

    use_coll = (switch_ps_time == 1 .or. (switch_ps_time == 4 .and. ps_tau_v > 0.0_r2))
    use_pr   = (switch_ps_time == 2 .or. switch_ps_time == 4)
    use_age  = (switch_ps_time == 3 .or. (switch_ps_time == 4 .and. ps_age_yr > 0.0_r2))

    if (allocated(ps_rloss)) deallocate(ps_rloss)
    allocate(ps_rloss(n_grain_radius, n_grain_components))
    ps_rloss = -1.0_r2
    if (allocated(ps_aloss)) deallocate(ps_aloss)
    allocate(ps_aloss(n_grain_radius, n_grain_components))
    ps_aloss = -1.0_r2
    s_h = 0.0_r2

    do ics = 1, n_grain_components
      if (ps_Y(ics) <= 0.0_r2) cycle
      ps_rloss(:, ics) = 0.0_r2
      icto  = ps_to_comp(ics)
      rho_w = ps_fice(ics) * particle_bulk_density(ics) * 1000.0_r2         ! [kg/m^3]
      ! distance entries: the parent orbit through the distance r. With orbit bins the distance
      ! entries are only used by circular belts, so the parent orbits there are circular.
      do idd = 1, n_dist
        r = dist_grid(idd)
        e_p = 0.0_r2;  a_p = r
        if (switch_ps_orbit == 2 .and. .not. orbit_mode) then
          e_p = ecc_parent_dist(idd)
          select case (switch_blowout_release)
          case (2);      a_p = r / (1.0_r2 - e_p)
          case (3);      a_p = r / (1.0_r2 + e_p)
          case default;  a_p = r;  e_p = 0.0_r2
          end select
        end if
        if (use_coll) s_h = blowout_radius(ics, idd)
        do igg = 1, n_grain_radius
          iss = sp_dist(igg, ics, idd)
          ! grains that are unbound when released leave within an orbit: blow-out removes them as
          ! they are, before they could lose their ice
          if (switch_blowout == 1 .and. beta_sp(igg, iss) >= beta_crit_dist(idd)) cycle
          if (ice_lost(igg, iss, ics, rho_w, r, a_p, e_p, r, s_h, Gq(igg, iss), Qm(igg, iss), &
                       use_coll, use_pr, use_age)) then
            ps_rloss(igg, ics) = max(ps_rloss(igg, ics), r)
            if (npm(igg, ics, idd) > 0.0_r2) then
              if (icto > 0) npm(igg, icto, idd) = npm(igg, icto, idd) + npm(igg, ics, idd)
              npm(igg, ics, idd) = 0.0_r2
            end if
          end if
        end do
      end do
      ! orbit bins of eccentric belts: parent orbit (a, e of the belt), released at the point given
      ! by blowout_release; Kepler frequency at a; species and s_1/2 at the release distance
      if (orbit_mode) then
        ! the distance entries serve the circular belts only: no ice-loss radius without them
        if (all(ob_count(1:n_belts) > 0)) ps_rloss(:, ics) = 0.0_r2
        ps_aloss(:, ics) = 0.0_r2
        do ob = 1, n_obins
          a_c = sqrt(ob_alo(ob) * ob_ahi(ob))
          e_p = belt_ecc(ob_belt(ob))
          select case (switch_blowout_release)
          case (2);      r = a_c * (1.0_r2 - e_p)
          case (3);      r = a_c * (1.0_r2 + e_p)
          case default;  r = a_c;  e_p = 0.0_r2
          end select
          idd = dist_index(r)
          if (use_coll) s_h = blowout_radius(ics, idd)
          do igg = 1, n_grain_radius
            iss = sp_dist(igg, ics, idd)
            if (ice_lost(igg, iss, ics, rho_w, r, a_c, e_p, a_c, s_h, Gq(igg, iss), Qm(igg, iss), &
                         use_coll, use_pr, use_age)) then
              ps_lost_ob(igg, ics, ob) = .true.
              ps_aloss(igg, ics) = max(ps_aloss(igg, ics), a_c)
            end if
          end do
        end do
      end if
      if (verbose) then
        nshow = max(1, n_grain_radius / 4)
        if (orbit_mode) then
          print *, "    component ", ics, ": ice lost by photosputtering on orbits out to (a [micron]: ", &
                   "semi-major axis [au]; circular belts: distance [au])"
          do ks = 1, n_grain_radius, nshow
            print *, "      ", real_to_string(rad(ks)), ": ", real_to_string(ps_aloss(ks, ics) / au), &
                     "; ", real_to_string(ps_rloss(ks, ics) / au)
          end do
        else
          print *, "    component ", ics, ": ice lost by photosputtering out to (a [micron]: r [au])"
          do ks = 1, n_grain_radius, nshow
            print *, "      ", real_to_string(rad(ks)), ": ", real_to_string(ps_rloss(ks, ics) / au)
          end do
          if (mod(n_grain_radius - 1, nshow) /= 0) print *, "      ", real_to_string(rad(n_grain_radius)), &
                  ": ", real_to_string(ps_rloss(n_grain_radius, ics) / au)
        end if
      end if
    end do

    if (verbose) then
      print *, "    photosputtering: band ", real_to_string(ps_lam1), " - ", real_to_string(ps_lam2), &
               " micron, ", real_to_string(ps_photons_1au * 1.0e-4_r2), " photons/cm^2/s at 1 au"
      if (ps_n_model > 0) print *, "    note: ps_uv_file does not cover ", ps_n_model, " of ", nb, &
                                   " band points; the stellar model is used there"
      if (ps_q_extrap .and. switch_ps_absorption == 1) &
        print *, "    note: the band extends beyond the optical data; Qabs at the nearest wavelength is used"
      if (n_uv == 0 .and. switch_BB_or_readingspec /= 2 .and. Tstar < 6000.0_r2) &
        print *, "    note: photospheric UV only; the chromospheric UV of cool stars is not included (ps_uv_file)"
    end if

    deallocate(Gq, Qm)
    if (allocated(uv_lam)) deallocate(uv_lam, uv_F)
    if (allocated(loglam)) deallocate(loglam)


  end subroutine photosputtering


  ! photosputtering: does a grain (size jg, species js) of component ics (ice mass per volume
  ! rho_w) lose its ice? Released at distance rr from a parent orbit (ap, ep) (local: rr = ap = r,
  ! ep = 0); Kepler frequency at rk; s_1/2 = sh; gq, qm: absorbed band photons at 1 au and band
  ! mean of Q of this grain
  logical function ice_lost(jg, js, ics, rho_w, rr, ap, ep, rk, sh, gq, qm, use_coll, use_pr, use_age)
    integer,  intent(in) :: jg, js, ics
    real(r2), intent(in) :: rho_w, rr, ap, ep, rk, sh, gq, qm
    logical,  intent(in) :: use_coll, use_pr, use_age
    real(r2), parameter  :: m_h2o = 18.01528_r2 * 1.66053906660e-27_r2       ! [kg]
    real(r2), parameter  :: year  = 365.25_r2 * 86400.0_r2                   ! [s]
    real(r2) :: b, inv_ag, a_g, ome2, f_orb, nabs, t_ps, t_ref, tc0, s, t_pr
    ice_lost = .false.
    b = beta_sp(jg, js)
    ! orbit of the grain (ps_orbit = 'release')
    f_orb = 1.0_r2
    a_g = rr;  ome2 = 1.0_r2
    if (switch_ps_orbit == 2) then
      inv_ag = (1.0_r2 / ap - 2.0_r2 * b / rr) / (1.0_r2 - b)
      if (b >= 1.0_r2 .or. inv_ag <= 0.0_r2) return              ! unbound: removed by blow-out
      a_g  = 1.0_r2 / inv_ag
      ome2 = ap * (1.0_r2 - ep**2) / ((1.0_r2 - b) * a_g)
      ome2 = min(max(ome2, 1.0e-12_r2), 1.0_r2)
      f_orb = rr**2 / (a_g**2 * sqrt(ome2))
    end if
    nabs = gq * (au / rr)**2 * f_orb + qm * ps_isrf_flux
    if (nabs <= 0.0_r2) return
    t_ps  = 4.0_r2 * rho_w * rad(jg) * micron / (m_h2o * ps_Y(ics) * nabs)
    t_ref = huge(1.0_r2)
    if (use_coll) then
      tc0 = 1.0_r2 / (ps_tau_v * sqrt(G_grav * Mstar / rk**3))
      s = rad(jg)
      if (s < 100.0_r2 * sh) then
        t_ref = min(t_ref, tc0 * ((s / (1.2_r2 * sh))**(-2) + (s / (100.0_r2 * sh))**2.7_r2))
      else
        t_ref = min(t_ref, tc0 * (s / (100.0_r2 * sh))**0.3_r2)
      end if
    end if
    if (use_pr .and. b > 0.0_r2) then
      if (switch_ps_orbit == 2) then
        t_pr = a_g**2 * c * ome2**1.5_r2 / (b * G_grav * Mstar * (2.0_r2 + 3.0_r2 * (1.0_r2 - ome2)))
      else
        t_pr = rr**2 * c / (2.0_r2 * b * G_grav * Mstar)
      end if
      t_ref = min(t_ref, t_pr)
    end if
    if (use_age) t_ref = min(t_ref, ps_age_yr * year)
    ice_lost = (t_ps < t_ref)
  end function ice_lost


  ! largest grain radius [micron] with beta >= 1/2 for component icomp at distance index idist
  ! (log interpolation between the size bins); the smallest grain if no grain reaches 1/2
  function blowout_radius(icomp, idist) result(s_half)
    integer, intent(in) :: icomp, idist
    real(r2)            :: s_half
    integer  :: j
    real(r2) :: b1, b2
    s_half = rad(1)
    do j = n_grain_radius, 1, -1
      b1 = beta_sp(j, sp_dist(j, icomp, idist))
      if (b1 >= 0.5_r2) then
        if (j == n_grain_radius) then
          s_half = rad(j)
        else
          b2 = beta_sp(j + 1, sp_dist(j + 1, icomp, idist))
          s_half = exp(log(rad(j)) + (b1 - 0.5_r2) / (b1 - b2) * (log(rad(j + 1)) - log(rad(j))))
        end if
        return
      end if
    end do
  end function blowout_radius


  pure function sort_key(icomp) result(key)
    integer, intent(in) :: icomp
    real(r2)            :: key
    key = T_subl(icomp)
    if (key <= 0.0_r2) key = huge(1.0_r2)
  end function sort_key

end subroutine particle_abundances
