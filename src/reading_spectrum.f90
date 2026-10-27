subroutine reading_spectrum

  !---------------------------------------------------------------------------------------------------
  ! Stellar radiation field
  !
  ! The star is described by its surface intensity I*(lambda) [W/m^2/m/sr]:
  !   black body : I* = B_lambda(T*)
  !   spectrum   : I* = F_lambda(d) d^2 / (pi R*^2),  F_lambda = F_nu c / lambda^2,  F_nu in Jy
  !
  ! Derived:
  !   L*                      = 4 pi^2 R*^2 int I* dlambda
  !   star_I_obs(il)          = I*(lambda_obs)
  !   stellar_flux_Jy(il)     = pi R*^2 I* / d^2 converted to Jy
  !
  ! Fixes with respect to the previous version: the spectrum arrays were deallocated right after
  ! reading, and Planck_fn returned 0 for T* = 0 even when a spectrum was read, so a measured
  ! spectrum produced no stellar heating at all. The radius of the star is now always required.
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2), allocatable, dimension(:,:) :: tab
  integer  :: i, il
  real(r2) :: integ

  if (switch_BB_or_readingspec >= 2) then

    if (switch_BB_or_readingspec == 2) then
      call read_table(trim(spec_star), 2, tab, n_spec)
      if (dms_ierr /= 0) return
      allocate(spec_lambda(n_spec), spec_I(n_spec))
      spec_lambda(:) = tab(1, :) * micron
      spec_I(:)      = tab(2, :) * 1.0e-26_r2 * c / spec_lambda(:)**2 &
                       * dist_star**2 / (pi * Rstar**2)
    else
      call read_stellar_grid
      if (dms_ierr /= 0) return
    end if

    do i = 2, n_spec
      if (spec_lambda(i) <= spec_lambda(i-1)) then
        call dms_fail("the wavelengths of the stellar spectrum must increase strictly")
        return
      end if
    end do

    ! luminosity: trapezoidal rule + analytic Rayleigh-Jeans tail
    integ = 0.0_r2
    do i = 1, n_spec - 1
      integ = integ + 0.5_r2 * (spec_I(i) + spec_I(i+1)) * (spec_lambda(i+1) - spec_lambda(i))
    end do
    integ = integ + spec_I(n_spec) * spec_lambda(n_spec) / 3.0_r2
    Lstar = 4.0_r2 * pi * pi * Rstar**2 * integ

    ! effective temperature equivalent (for information only; the grid keeps the requested T_eff)
    if (switch_BB_or_readingspec == 2) then
      Tstar = (Lstar / (4.0_r2 * pi * Rstar**2 * sigma_SB))**0.25_r2
    else if (verbose) then
      print *, "    T_eff from the integrated model spectrum: ", &
               real_to_string((Lstar / (4.0_r2 * pi * Rstar**2 * sigma_SB))**0.25_r2), " K"
    end if

    if (verbose) print *, " "
    if (verbose) print *, "    stellar spectrum: ", n_spec, " points, ", &
             real_to_string(spec_lambda(1)/micron), " - ", &
             real_to_string(spec_lambda(n_spec)/micron), " micron"

  else

    Lstar = 4.0_r2 * pi * Rstar**2 * sigma_SB * Tstar**4

  end if

  if (verbose) print *, "    L* = ", real_to_scientific_string(Lstar / 3.828e26_r2), " L_sun,  T_eff = ", &
           real_to_string(Tstar), " K"

  allocate(star_I_obs(n_wavelength_obs), stellar_flux_Jy(n_wavelength_obs))
  do il = 1, n_wavelength_obs
    star_I_obs(il)      = stellar_radiance(lambda_obs_m(il))
    stellar_flux_Jy(il) = pi * Rstar**2 * star_I_obs(il) / dist_star**2 &
                          * lambda_obs_m(il)**2 / c * 1.0e26_r2
  end do

contains

  !-------------------------------------------------------------------------------------------------
  ! Model-atmosphere grid (PHOENIX, ATLAS9/Castelli-Kurucz, ... prepared with
  ! python -m dms.stellar or by hand). The directory contains grid_index.txt with lines
  !     T_eff [K]   log g [cgs]   file name
  ! and one file per model with two columns: wavelength [micron], surface flux F_lambda
  ! [W m^-2 m^-1] (= pi I). The spectrum is interpolated bilinearly in (T_eff, log g) in log F_lambda,
  ! on the wavelength grid of the first model used.
  !-------------------------------------------------------------------------------------------------
  subroutine read_stellar_grid

    integer, parameter :: cap = 50000
    real(r2), allocatable :: tg(:), gg(:), tab1(:,:), lflux(:), lam0(:)
    character(len=512), allocatable :: fg(:)
    character(len=1024) :: line, iname
    integer  :: u, ios, n, k, it, ig2, j, nrow, nsel
    real(r2) :: tt, gv, t_lo, t_hi, wt, g_lo, g_hi, wg, w
    real(r2) :: tnode(2), wnode(2)
    character(len=512) :: fn
    logical  :: ex

    iname = trim(stellar_grid_dir)//"/grid_index.txt"
    inquire(file=trim(iname), exist=ex)
    if (.not. ex) then
      call dms_fail("model-atmosphere grid index not found: "//trim(iname));  return
    end if

    allocate(tg(cap), gg(cap), fg(cap))
    n = 0
    open(newunit=u, file=trim(iname), action="read", status="old", form="formatted")
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read(line, *, iostat=ios) tt, gv, fn
      if (ios /= 0) cycle
      if (n == cap) exit
      n = n + 1
      tg(n) = tt;  gg(n) = gv;  fg(n) = fn
    end do
    close(u)
    if (n == 0) then
      call dms_fail("no models listed in "//trim(iname));  return
    end if

    ! bracketing temperatures
    t_lo = -huge(1.0_r2);  t_hi = huge(1.0_r2)
    do k = 1, n
      if (tg(k) <= Tstar .and. tg(k) > t_lo) t_lo = tg(k)
      if (tg(k) >= Tstar .and. tg(k) < t_hi) t_hi = tg(k)
    end do
    if (t_lo < 0.0_r2 .or. t_hi > 1.0e29_r2) then
      call dms_fail("T_eff = "//real_to_string(Tstar)//" K is outside the model grid");  return
    end if
    tnode = [t_lo, t_hi]
    wt = 0.0_r2
    if (t_hi > t_lo) wt = (Tstar - t_lo) / (t_hi - t_lo)
    wnode = [1.0_r2 - wt, wt]

    nsel = 0
    do it = 1, 2
      if (wnode(it) <= 0.0_r2) cycle
      ! bracketing log g at this temperature
      g_lo = -huge(1.0_r2);  g_hi = huge(1.0_r2)
      do k = 1, n
        if (abs(tg(k) - tnode(it)) > 1.0e-6_r2) cycle
        if (gg(k) <= logg_star .and. gg(k) > g_lo) g_lo = gg(k)
        if (gg(k) >= logg_star .and. gg(k) < g_hi) g_hi = gg(k)
      end do
      if (g_lo < -1.0e29_r2 .or. g_hi > 1.0e29_r2) then
        call dms_fail("log g = "//real_to_string(logg_star)//" is outside the model grid at T_eff = "// &
                      real_to_string(tnode(it))//" K");  return
      end if
      wg = 0.0_r2
      if (g_hi > g_lo) wg = (logg_star - g_lo) / (g_hi - g_lo)
      do ig2 = 1, 2
        if (ig2 == 1) then
          gv = g_lo;  w = wnode(it) * (1.0_r2 - wg)
        else
          gv = g_hi;  w = wnode(it) * wg
        end if
        if (w <= 0.0_r2) cycle
        do k = 1, n
          if (abs(tg(k) - tnode(it)) <= 1.0e-6_r2 .and. abs(gg(k) - gv) <= 1.0e-6_r2) exit
        end do
        call read_table(trim(stellar_grid_dir)//"/"//trim(fg(k)), 2, tab1, nrow)
        if (dms_ierr /= 0) return
        if (nsel == 0) then
          allocate(lam0(nrow), lflux(nrow))
          lam0  = tab1(1, :)
          lflux = 0.0_r2
        end if
        do j = 1, size(lam0)
          lflux(j) = lflux(j) + w * interp_lin(log(tab1(1, :)), log(max(tab1(2, :), 1.0e-300_r2)), log(lam0(j)))
        end do
        nsel = nsel + 1
        if (verbose) print *, "    model atmosphere: ", trim(fg(k)), "  weight ", real_to_string(w)
      end do
    end do

    n_spec = size(lam0)
    allocate(spec_lambda(n_spec), spec_I(n_spec))
    spec_lambda(:) = lam0(:) * micron
    spec_I(:)      = exp(lflux(:)) / pi
    where (spec_I < 1.0e-250_r2) spec_I = 0.0_r2

  end subroutine read_stellar_grid

end subroutine reading_spectrum
