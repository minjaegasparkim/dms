subroutine temperature_distribution_cals

  !---------------------------------------------------------------------------------------------------
  ! Equilibrium temperature of every grain (size ig, component ic) at every distance of the
  ! distance grid (n_dist points between the pericentre Rmin(1-e) and the apocentre Rmax(1+e)).
  !
  ! Energy balance of an optically thin grain (Wolf & Hillenbrand 2003, eq. 5):
  !     r(T) = R*/2 * sqrt( int Qabs I*(lambda) dlambda / int Qabs B_lambda(T) dlambda )
  !
  ! Method (replaces the former iterative search of T_max/T_min in 0.5 K steps and the linear
  ! search through the temperature table for every grain and every volume element):
  !   1. tabulate r(T) for n_temp_table temperatures between 1 K and 5000 K (log-spaced), for all
  !      grains at once (matrix product Qabs^T x B(T)),
  !   2. invert r(T) -> T(r) on the distance grid by interpolation in log r / log T,
  !   3. look-up during the integration: O(1) via the index of the distance grid.
  ! With temperature-dependent optical constants, the data set with the temperature tag closest
  ! to T is used at each T of the table (r(T) is forced to be monotonic afterwards).
  !
  ! Also computed here: <Qabs>* = int Qabs I* / int I* (for the fractional luminosity).
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2), parameter :: T_table_min = 1.0_r2, T_table_max = 5000.0_r2

  real(r2), allocatable, dimension(:)     :: T_table, Istar_ri
  real(r2), allocatable, dimension(:,:)   :: Bw, Sabs, r_tab
  real(r2), allocatable, dimension(:,:,:) :: Dall
  integer,  allocatable, dimension(:)     :: sp_of_T
  integer  :: k, ilr, ig, ic, is, id, i, nrow, step, unit
  real(r2) :: lr, t, Tg
  real(r2), allocatable, dimension(:,:) :: logr

  if (verbose) print *, " "
  if (verbose) print *, " 5) Radial temperature distribution"
  if (verbose) print *, " "

  allocate(T_table(n_temp_table), Istar_ri(n_wavelength_int))
  allocate(Bw(n_wavelength_int, n_temp_table))
  allocate(Sabs(n_grain_radius, n_species), Dall(n_grain_radius, n_temp_table, n_species))
  allocate(Qstar_mean(n_grain_radius, n_species))

  do k = 1, n_temp_table
    T_table(k) = T_table_min * (T_table_max / T_table_min)**(real(k-1, kind=r2) / real(n_temp_table-1, kind=r2))
  end do

  do ilr = 1, n_wavelength_int
    Istar_ri(ilr) = stellar_radiance(wavelength_int(ilr) * micron)
  end do

  ! B_lambda(T) times the integration weights
  !$omp parallel do private(k, ilr)
  do k = 1, n_temp_table
    do ilr = 1, n_wavelength_int
      Bw(ilr, k) = Planck_fn(T_table(k), wavelength_int(ilr) * micron) * dlambda_int(ilr)
    end do
  end do
  !$omp end parallel do

  do is = 1, n_species
    do ig = 1, n_grain_radius
      Sabs(ig, is)       = sum(Qabs_ri(:, ig, is) * Istar_ri(:) * dlambda_int(:))
      Qstar_mean(ig, is) = Sabs(ig, is) / max(sum(Istar_ri(:) * dlambda_int(:)), tiny(1.0_r2))
    end do
    Dall(:, :, is) = matmul(transpose(Qabs_ri(:, :, is)), Bw)      ! (ig, k)
  end do

  !------------------------------------------------------------------------------------------------
  ! r(T) for every grain and component (species chosen by temperature)
  !------------------------------------------------------------------------------------------------

  dist_lo = r_in_eff
  dist_hi = r_max_3d
  if (dist_hi <= dist_lo) dist_hi = dist_lo * (1.0_r2 + 1.0e-6_r2)
  ! logarithmic grid: constant relative resolution (also for disks reaching close to the star)
  dist_step = log(dist_hi / dist_lo) / real(n_dist - 1, kind=r2)

  allocate(dist_grid(n_dist))
  do id = 1, n_dist
    dist_grid(id) = dist_lo * exp(real(id - 1, kind=r2) * dist_step)
  end do

  allocate(T_dist(n_grain_radius, n_grain_components, n_dist), &
           sp_dist(n_grain_radius, n_grain_components, n_dist), T1_dist(n_dist))
  allocate(r_tab(n_temp_table, n_grain_radius), sp_of_T(n_temp_table), logr(n_temp_table, n_grain_radius))

  do ic = 1, n_grain_components

    ! species used at each temperature of the table
    do k = 1, n_temp_table
      sp_of_T(k) = nearest_species(ic, T_table(k))
    end do

    do ig = 1, n_grain_radius
      do k = 1, n_temp_table
        is = sp_of_T(k)
        if (Dall(ig, k, is) > 0.0_r2 .and. Sabs(ig, is) > 0.0_r2) then
          r_tab(k, ig) = 0.5_r2 * Rstar * sqrt(Sabs(ig, is) / Dall(ig, k, is))
        else
          r_tab(k, ig) = huge(1.0_r2) / 1.0e10_r2
        end if
      end do
      ! r must decrease with increasing T
      do k = n_temp_table - 1, 1, -1
        r_tab(k, ig) = max(r_tab(k, ig), r_tab(k+1, ig) * (1.0_r2 + 1.0e-12_r2))
      end do
      logr(:, ig) = log(r_tab(:, ig))
    end do

    ! invert: T(r) on the distance grid
    !$omp parallel do private(id, ig, lr, i, t, Tg)
    do id = 1, n_dist
      lr = log(dist_grid(id))
      do ig = 1, n_grain_radius
        if (lr >= logr(1, ig)) then
          Tg = T_table(1)
        else if (lr <= logr(n_temp_table, ig)) then
          Tg = T_table(n_temp_table)
        else
          i  = locate(logr(:, ig), lr)
          t  = (lr - logr(i, ig)) / (logr(i+1, ig) - logr(i, ig))
          Tg = exp(log(T_table(i)) + t * (log(T_table(i+1)) - log(T_table(i))))
        end if
        T_dist(ig, ic, id)  = Tg
        sp_dist(ig, ic, id) = nearest_species(ic, Tg)
      end do
    end do
    !$omp end parallel do

  end do

  T1_dist(:) = T_dist(1, 1, :)

  if (verbose) print *, "    distance grid: ", n_dist, " points, ", real_to_string(dist_lo/au), " - ", &
           real_to_string(dist_hi/au), " au"
  if (verbose) print *, "    smallest grains, component 1: T = ", real_to_string(T_dist(1,1,1)), " K at ", &
           real_to_string(dist_lo/au), " au,  ", real_to_string(T_dist(1,1,n_dist)), " K at ", &
           real_to_string(dist_hi/au), " au"
  if (verbose) print *, "    largest  grains, component 1: T = ", real_to_string(T_dist(n_grain_radius,1,1)), &
           " K at ", real_to_string(dist_lo/au), " au,  ", &
           real_to_string(T_dist(n_grain_radius,1,n_dist)), " K at ", real_to_string(dist_hi/au), " au"

  !------------------------------------------------------------------------------------------------
  ! write the temperature distribution (at most ~2000 rows)
  !------------------------------------------------------------------------------------------------

  if (switch_temperature_distribution == 1 .and. write_output) then

    call make_dir(results_dir)
    open(newunit=unit, file=trim(results_dir)//"Temperature_Distribution"//file_suffix(0)//".dat", &
         action="write", status="replace", form="formatted")
    write(unit, '(a)') "# distance [au], then T [K] for each component (outer loop) and grain radius"
    write(unit, '(a)', advance='no') "# radii [micron]:"
    do ig = 1, n_grain_radius
      write(unit, '(es11.3)', advance='no') rad(ig)
    end do
    write(unit, *)
    step = max(1, n_dist / 2000)
    nrow = 0
    do id = 1, n_dist, step
      write(unit, '(es14.6)', advance='no') dist_grid(id) / au
      do ic = 1, n_grain_components
        do ig = 1, n_grain_radius
          write(unit, '(es12.4)', advance='no') T_dist(ig, ic, id)
        end do
      end do
      write(unit, *)
      nrow = nrow + 1
    end do
    close(unit)
    if (verbose) print *, "    temperature distribution written (", nrow, " distances)"

  end if

  deallocate(T_table, Istar_ri, Bw, Sabs, Dall, r_tab, sp_of_T, logr)

contains

  pure function nearest_species(icomp, temp) result(isp)
    integer, intent(in)  :: icomp
    real(r2), intent(in) :: temp
    integer              :: isp, kk
    real(r2)             :: dbest
    isp = comp_first_species(icomp)
    if (comp_n_species(icomp) == 1) return
    dbest = huge(1.0_r2)
    do kk = comp_first_species(icomp), comp_first_species(icomp) + comp_n_species(icomp) - 1
      if (abs(species_T(kk) - temp) < dbest) then
        dbest = abs(species_T(kk) - temp)
        isp   = kk
      end if
    end do
  end function nearest_species

end subroutine temperature_distribution_cals
