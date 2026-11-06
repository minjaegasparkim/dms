subroutine reading_temperature_dependent_data

  !---------------------------------------------------------------------------------------------------
  ! Reads the optical constants of all dust components and sets up the "optical species".
  !
  ! A component without temperature-dependent optical constants has one data set (its optical-data
  ! file) and one species. A component with temperature-dependent data (DMS 2.1: tdep = .true.,
  ! any number of components) has the data sets <base>_<T>K[.ext] found in the optical-data
  ! directory (module tdep_search), or the sets listed with tdep_temperatures/tdep_files (input of
  ! DMS 2.0, one component). The species of such a component are
  !   tdep_mode = 'nearest' : one species per data set; a grain uses the set whose temperature is
  !                           closest to its own equilibrium temperature;
  !   tdep_mode = 'linear'  : in addition, species between neighbouring sets with n and k
  !                           interpolated linearly in T, spaced by at most tdep_dt; a grain uses
  !                           the closest of these, i.e. n, k interpolated to within tdep_dt/2 of
  !                           its temperature.
  ! The choice is self-consistent (also in the temperature calculation) and is made for every grain
  ! size, component and distance (temperature_distribution_cals).
  !
  ! The wavelength grid of the optical data ("ri grid") is the union of the wavelengths of all
  ! files, restricted to the range covered by all files; every data set is interpolated onto it
  ! (linear in log lambda). (Before this fix it was the grid of the first file only, so that a
  ! finely sampled band of another component, e.g. the 3 micron ice band in an external mixture
  ! listed after a coarsely sampled refractory material, was smoothed.)
  ! The files may have any number of lines; headers and comment lines are skipped.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions
  use tdep_search

  implicit none

  type table_t
    real(r2), allocatable, dimension(:,:) :: t
    integer :: n
  end type table_t

  type(table_t), allocatable, dimension(:) :: tabs
  integer  :: ic, is, k, kk, ia, ilr, n0, i1, i2, n_sets, iset, j, m, nf
  real(r2) :: lam_lo, lam_hi, lam_new, dT
  real(r2), allocatable, dimension(:) :: grid_all, grid_new, tf
  real(r2), allocatable, dimension(:,:) :: n_set, k_set
  character(len=512), allocatable, dimension(:) :: ff
  character(len=512) :: errmsg
  logical :: is_tdep

  if (verbose) print *, " "
  if (verbose) print *, " 2) Reading optical data"
  if (verbose) print *, " "

  !------------------------------------------------------------------------------------------------
  ! data sets of every component
  !------------------------------------------------------------------------------------------------

  if (.not. allocated(comp_tdep)) then
    allocate(comp_tdep(n_grain_components))
    comp_tdep(:) = .false.
  end if
  if (switch_Tdep_optics == 1) comp_tdep(Tdep_comp) = .true.        ! input of DMS 2.0

  allocate(comp_first_set(n_grain_components), comp_n_sets(n_grain_components))
  allocate(set_T(0), set_file(0))
  n_sets = 0
  do ic = 1, n_grain_components
    comp_first_set(ic) = n_sets + 1
    if (.not. comp_tdep(ic)) then
      call add_set(-1.0_r2, dustfile_name(ic))
    else if (switch_Tdep_optics == 1 .and. ic == Tdep_comp) then
      ! DMS 2.0: temperatures and files given explicitly (sorted by temperature)
      allocate(tf(n_Tdep), ff(n_Tdep))
      tf(:) = Tdep_T(1:n_Tdep);  ff(:) = Tdep_file(1:n_Tdep)
      call sort_sets(tf, ff)
      do k = 2, n_Tdep
        if (abs(tf(k) - tf(k-1)) <= 1.0e-9_r2 * tf(k)) then
          call dms_fail("tdep_temperatures: the temperature "//real_to_string(tf(k))//" K is given twice")
          return
        end if
      end do
      do k = 1, n_Tdep
        call add_set(tf(k), ff(k))
      end do
      deallocate(tf, ff)
    else
      ! DMS 2.1: files <base>_<T>K[.ext] in the optical-data directory
      call find_tdep_files(dustfile_name(ic), ri_dir, nf, tf, ff, errmsg)
      if (len_trim(errmsg) > 0) then
        call dms_fail(trim(errmsg));  return
      end if
      if (nf == 0) then
        call dms_fail("component "//integer_to_string(ic)//": no files "//trim(dustfile_name(ic))// &
                      "_<T>K[.ext] in "//trim(ri_dir)//" (with tdep = .true. optical_data is the base "// &
                      "name without _<T>K)")
        return
      else if (nf == 1) then
        if (verbose) print *, "    warning: component ", ic, ": only one file ", trim(ff(1)), &
                              " found; it is used without temperature dependence"
        comp_tdep(ic) = .false.
        call add_set(-1.0_r2, ff(1))
      else
        do k = 1, nf
          call add_set(tf(k), ff(k))
        end do
      end if
      deallocate(tf, ff)
    end if
    comp_n_sets(ic) = n_sets - comp_first_set(ic) + 1
  end do

  !------------------------------------------------------------------------------------------------
  ! read all files
  !------------------------------------------------------------------------------------------------

  allocate(tabs(n_sets))
  do iset = 1, n_sets
    call read_table(trim(ri_dir)//trim(set_file(iset)), 3, tabs(iset)%t, tabs(iset)%n)
    if (dms_ierr /= 0) return
    if (tabs(iset)%t(1, tabs(iset)%n) < tabs(iset)%t(1, 1)) then          ! make ascending
      tabs(iset)%t = tabs(iset)%t(:, tabs(iset)%n:1:-1)
    end if
  end do

  lam_lo = -huge(1.0_r2)
  lam_hi =  huge(1.0_r2)
  do iset = 1, n_sets
    lam_lo = max(lam_lo, tabs(iset)%t(1, 1))
    lam_hi = min(lam_hi, tabs(iset)%t(1, tabs(iset)%n))
  end do

  if (lam_hi <= lam_lo) then
    call dms_fail("the optical-data files have no common wavelength range")
    return
  end if

  ! ri grid = union of the wavelengths of all files within [lam_lo, lam_hi]
  n0 = 0
  do iset = 1, n_sets
    n0 = n0 + tabs(iset)%n
  end do
  allocate(grid_all(n0), grid_new(n0))
  n_wavelength_int = 0
  do iset = 1, n_sets
    is = iset
    ! wavelengths of this file inside the common range (ascending)
    i1 = 1
    do while (i1 <= tabs(is)%n)
      if (tabs(is)%t(1, i1) >= lam_lo * (1.0_r2 - 1.0e-12_r2)) exit
      i1 = i1 + 1
    end do
    i2 = tabs(is)%n
    do while (i2 >= 1)
      if (tabs(is)%t(1, i2) <= lam_hi * (1.0_r2 + 1.0e-12_r2)) exit
      i2 = i2 - 1
    end do
    if (i2 < i1) cycle
    ! merge with the grid so far, dropping duplicates (relative difference < 1e-9)
    ia = 1;  k = i1;  kk = 0
    do while (ia <= n_wavelength_int .or. k <= i2)
      if (k > i2) then
        lam_new = grid_all(ia);  ia = ia + 1
      else if (ia > n_wavelength_int) then
        lam_new = tabs(is)%t(1, k);  k = k + 1
      else if (grid_all(ia) <= tabs(is)%t(1, k)) then
        lam_new = grid_all(ia);  ia = ia + 1
      else
        lam_new = tabs(is)%t(1, k);  k = k + 1
      end if
      if (kk > 0) then
        if (abs(lam_new - grid_new(kk)) <= 1.0e-9_r2 * lam_new) cycle
      end if
      kk = kk + 1
      grid_new(kk) = lam_new
    end do
    n_wavelength_int = kk
    grid_all(1:kk) = grid_new(1:kk)
  end do

  allocate(wavelength_int(n_wavelength_int), dlambda_int(n_wavelength_int))
  wavelength_int(:) = grid_all(1:n_wavelength_int)
  deallocate(grid_all, grid_new)

  allocate(n_set(n_wavelength_int, n_sets), k_set(n_wavelength_int, n_sets))
  do iset = 1, n_sets
    do ilr = 1, n_wavelength_int
      n_set(ilr, iset) = interp_lin(log(tabs(iset)%t(1,:)), tabs(iset)%t(2,:), log(wavelength_int(ilr)))
      k_set(ilr, iset) = interp_lin(log(tabs(iset)%t(1,:)), tabs(iset)%t(3,:), log(wavelength_int(ilr)))
    end do
  end do

  !------------------------------------------------------------------------------------------------
  ! species
  !------------------------------------------------------------------------------------------------

  allocate(comp_first_species(n_grain_components), comp_n_species(n_grain_components))
  n_species = 0
  do ic = 1, n_grain_components
    comp_first_species(ic) = n_species + 1
    if (comp_n_sets(ic) > 1 .and. switch_tdep_mode == 2) then
      do k = 1, comp_n_sets(ic) - 1
        iset = comp_first_set(ic) + k - 1
        n_species = n_species + max(1, ceiling((set_T(iset+1) - set_T(iset)) / tdep_step - 1.0e-9_r2))
      end do
      n_species = n_species + 1
    else
      n_species = n_species + comp_n_sets(ic)
    end if
    comp_n_species(ic) = n_species - comp_first_species(ic) + 1
  end do

  allocate(species_comp(n_species), species_T(n_species), species_file(n_species), &
           species_is_data(n_species), species_set_lo(n_species), species_set_hi(n_species), &
           species_w(n_species))
  allocate(n_ri(n_wavelength_int, n_species), k_ri(n_wavelength_int, n_species))

  is = 0
  do ic = 1, n_grain_components
    is_tdep = comp_n_sets(ic) > 1
    do k = 1, comp_n_sets(ic)
      iset = comp_first_set(ic) + k - 1
      m = 1
      if (is_tdep .and. switch_tdep_mode == 2 .and. k < comp_n_sets(ic)) &
        m = max(1, ceiling((set_T(iset+1) - set_T(iset)) / tdep_step - 1.0e-9_r2))
      do j = 0, m - 1
        is = is + 1
        species_comp(is) = ic
        species_set_lo(is) = iset
        if (j == 0) then
          species_set_hi(is) = iset;  species_w(is) = 0.0_r2;  species_is_data(is) = .true.
          species_T(is) = set_T(iset);  species_file(is) = set_file(iset)
          n_ri(:, is) = n_set(:, iset);  k_ri(:, is) = k_set(:, iset)
        else
          species_set_hi(is) = iset + 1;  species_w(is) = real(j, r2) / real(m, r2)
          species_is_data(is) = .false.
          dT = set_T(iset+1) - set_T(iset)
          species_T(is) = set_T(iset) + species_w(is) * dT
          species_file(is) = "(interpolated)"
          n_ri(:, is) = (1.0_r2 - species_w(is)) * n_set(:, iset) + species_w(is) * n_set(:, iset+1)
          k_ri(:, is) = (1.0_r2 - species_w(is)) * k_set(:, iset) + species_w(is) * k_set(:, iset+1)
        end if
      end do
    end do
    if (.not. is_tdep) species_T(comp_first_species(ic)) = -1.0_r2
  end do

  if (verbose) then
    do ic = 1, n_grain_components
      if (comp_n_sets(ic) == 1) then
        print *, "    component ", ic, ": ", trim(set_file(comp_first_set(ic))), "  (", &
                 tabs(comp_first_set(ic))%n, " wavelengths)"
      else
        print *, "    component ", ic, ": temperature-dependent data, ", comp_n_sets(ic), " sets, ", &
                 real_to_string(set_T(comp_first_set(ic))), " - ", &
                 real_to_string(set_T(comp_first_set(ic) + comp_n_sets(ic) - 1)), " K"
        do k = 1, comp_n_sets(ic)
          iset = comp_first_set(ic) + k - 1
          print *, "        ", real_to_string(set_T(iset)), " K  ", trim(set_file(iset)), "  (", &
                   tabs(iset)%n, " wavelengths)"
        end do
        if (switch_tdep_mode == 2) then
          print *, "        'linear': ", comp_n_species(ic), " species (n, k interpolated in T, step <= ", &
                   real_to_string(tdep_step), " K)"
        else
          print *, "        'nearest': each grain uses the set closest to its temperature"
        end if
      end if
    end do
    if (n_species > 300) print *, "    note: ", n_species, " optical species; the Mie calculation takes ", &
                                  "correspondingly longer (save_optics = .true. keeps them)"
  end if

  deallocate(n_set, k_set)

  ! trapezoidal integration weights [m]
  do ilr = 1, n_wavelength_int
    if (n_wavelength_int == 1) then
      dlambda_int(ilr) = 0.0_r2
    else if (ilr == 1) then
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(2) - wavelength_int(1)) * micron
    else if (ilr == n_wavelength_int) then
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(ilr) - wavelength_int(ilr-1)) * micron
    else
      dlambda_int(ilr) = 0.5_r2 * (wavelength_int(ilr+1) - wavelength_int(ilr-1)) * micron
    end if
  end do

  if (verbose) print *, " "
  if (verbose) print *, "    common wavelength grid: ", n_wavelength_int, " points, ", &
           real_to_string(wavelength_int(1)), " - ", real_to_string(wavelength_int(n_wavelength_int)), &
           " micron"

  !------------------------------------------------------------------------------------------------
  ! checks
  !------------------------------------------------------------------------------------------------

  if (minval(wavelength_obs) < wavelength_int(1) * (1.0_r2 - 1.0e-9_r2) .or. &
      maxval(wavelength_obs) > wavelength_int(n_wavelength_int) * (1.0_r2 + 1.0e-9_r2)) then
    call dms_fail("the observing wavelengths must lie within the wavelength range of the optical data ("// &
                  real_to_string(wavelength_int(1))//" - "//real_to_string(wavelength_int(n_wavelength_int))// &
                  " micron)")
    return
  end if

  if (stellar_radiance(wavelength_int(1) * micron) > 0.0_r2 .and. &
      wavelength_int(1) > 0.3_r2 * 2897.77_r2 / max(Tstar, 1.0_r2)) then
    if (verbose) print *, "    warning: the optical data start at ", real_to_string(wavelength_int(1)), &
             " micron; part of the stellar UV/optical flux is not absorbed by the grains."
  end if

  deallocate(tabs)

contains

  subroutine add_set(temp, fname)
    real(r2), intent(in)         :: temp
    character(len=*), intent(in) :: fname
    real(r2), allocatable, dimension(:) :: t2
    character(len=512), allocatable, dimension(:) :: f2
    allocate(t2(n_sets + 1), f2(n_sets + 1))
    if (n_sets > 0) then
      t2(1:n_sets) = set_T(1:n_sets);  f2(1:n_sets) = set_file(1:n_sets)
    end if
    n_sets = n_sets + 1
    t2(n_sets) = temp;  f2(n_sets) = fname
    call move_alloc(t2, set_T);  call move_alloc(f2, set_file)
  end subroutine add_set

  subroutine sort_sets(t, f)
    real(r2), dimension(:), intent(inout)         :: t
    character(len=512), dimension(:), intent(inout) :: f
    integer  :: i, jj
    real(r2) :: ts
    character(len=512) :: fs
    do i = 2, size(t)
      ts = t(i);  fs = f(i);  jj = i - 1
      do while (jj >= 1)
        if (t(jj) <= ts) exit
        t(jj+1) = t(jj);  f(jj+1) = f(jj);  jj = jj - 1
      end do
      t(jj+1) = ts;  f(jj+1) = fs
    end do
  end subroutine sort_sets

end subroutine reading_temperature_dependent_data
