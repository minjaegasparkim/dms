!__________________________________________________________________________________________|
!                                                                                          |
! 4th module : functions used in the whole code                                            |
!__________________________________________________________________________________________|
!==========================================================================================|
!  - coordinate transformations  detector <-> sky <-> disk frame                           |
!  - analytical density distributions (incl. eccentric discs)                              |
!  - Planck function and stellar surface intensity                                         |
!  - Henyey-Greenstein phase functions                                                     |
!  - distance-grid look-up, interpolation, robust table reader                             |
!  - file-name helpers                                                                     |
!                                                                                          |
! All functions here only READ global variables. They are thread-safe.                     |
!==========================================================================================|

module functions

  use datatype
  use constants
  use parameters
  use user_density

  implicit none

contains

  !_______________________________________________________________________________________|
  ! 0. Error handling                                                                      |
  !    Command-line program: print and stop. Python interface (library_mode): set the      |
  !    error flag; the calling routine returns and the driver stops the run.              |
  !_______________________________________________________________________________________|

  subroutine dms_fail(msg)

    character(len=*), intent(in) :: msg

    if (dms_ierr == 0) then
      dms_ierr   = 1
      dms_errmsg = msg
    end if
    if (verbose .or. .not. library_mode) write(*, '(a)') " <!> ERROR: "//trim(msg)
    if (.not. library_mode) error stop 1

  end subroutine dms_fail


  subroutine make_dir(path)                                 ! create a directory (with parents)

    character(len=*), intent(in) :: path
    logical :: ex

    inquire(file=trim(path)//"/.", exist=ex)
    if (.not. ex) call execute_command_line("mkdir -p '"//trim(path)//"'")

  end subroutine make_dir

  !_______________________________________________________________________________________|
  ! 1. Coordinate transformations                                                          |
  !_______________________________________________________________________________________|
  ! detector (xd, yd, z) -> sky (x, y, z):  rotation by the position angle about z         |
  !   x =  xd cos(PA) + yd sin(PA),   y = -xd sin(PA) + yd cos(PA)                         |
  ! sky -> disk (X, Y, Z): rotation by the inclination about the y axis                    |
  !   X = x cos(i) - z sin(i),   Y = y,   Z = x sin(i) + z cos(i)                          |
  ! For PA = 0 this is identical to the transformation of the previous DMS version.        |
  !_______________________________________________________________________________________|

  pure subroutine detector_to_disk(xd, yd, z, X, Y, Zd)

    real(r2), intent(in)  :: xd, yd, z
    real(r2), intent(out) :: X, Y, Zd
    real(r2)              :: xs, ys

    xs =  xd * cos_pa + yd * sin_pa
    ys = -xd * sin_pa + yd * cos_pa

    X  = xs * cos_inc - z * sin_inc
    Y  = ys
    Zd = xs * sin_inc + z * cos_inc

  end subroutine detector_to_disk


  pure subroutine disk_to_detector(X, Y, Zd, xd, yd, z)

    real(r2), intent(in)  :: X, Y, Zd
    real(r2), intent(out) :: xd, yd, z
    real(r2)              :: xs, ys

    xs =  X * cos_inc + Zd * sin_inc
    ys =  Y
    z  = -X * sin_inc + Zd * cos_inc

    xd = xs * cos_pa - ys * sin_pa
    yd = xs * sin_pa + ys * cos_pa

  end subroutine disk_to_detector

  !_______________________________________________________________________________________|
  ! 2. Density distributions                                                               |
  !_______________________________________________________________________________________|
  ! The disk consists of one or more belts (n_belts). Each belt k has a radial law          |
  ! (belt_law), inner/outer semi-major axes, eccentricity e and argument of pericentre      |
  ! omega (angle in the disk plane from the line of nodes +Y towards +X).                   |
  !                                                                                        |
  ! Eccentric belts: confocal orbits. A point at distance s and true anomaly f lies on the  |
  ! orbit with semi-major axis  a = s (1 + e cos f) / (1 - e^2).  For such a set of orbits  |
  ! the surface density is Sigma(s, f) ~ Sigma_0(a): the time a particle spends per unit   |
  ! azimuth (~ s^2) cancels the Jacobian s ds/da (= s^2/a), so the density becomes          |
  !        rho(s, f) = rho_0(a) * a / s,                                                   |
  ! which reduces to rho_0(s) for e = 0.                                                    |
  !                                                                                        |
  ! Vertical structure (common to all belts):                                              |
  !   wedge    : rho = rho_0 within |Z| <= sin(theta) r, s = r (3-d distance)              |
  !   Gaussian : rho = rho_0 exp(-Z^2 / 2h^2), h = (h/r) varpi, cut at |Z| = gauss_cut h,  |
  !              s = varpi (distance in the disk plane); rho_0 is the mid-plane density    |
  ! Law 5 calls user_density_fn (user_density.f90) instead of a built-in radial law.        |
  !_______________________________________________________________________________________|

  pure function radial_profile_fn(a, k) result(prof)

    real(r2), intent(in) :: a                                ! semi-major axis [m]
    integer,  intent(in) :: k                                ! belt
    real(r2)             :: prof
    real(r2)             :: s

    select case (belt_law(k))

    case (1)                                                 ! n ~ a^(-alpha)
      prof = (a / au)**(-belt_alpha(k))

    case (2)                                                 ! Ertel et al. (2012)
      if (a < belt_Rstart(k)) then
        prof = 0.0_r2
      else
        s    = (a - belt_Rstart(k)) / (belt_Rpeak(k) - belt_Rstart(k))
        prof = s**belt_alpha(k) * exp(belt_alpha(k) * (1.0_r2 - s))
      end if

    case (3)                                                 ! a^alpha * Gaussian
      prof = (a / au)**belt_alpha(k) * exp(-0.5_r2 * ((a - belt_mu(k)) / belt_sigma(k))**2)

    case (4)                                                 ! smooth broken power law
      s    = a / belt_Rc(k)                                  ! (Augereau et al. 1999)
      prof = 1.0_r2 / sqrt(s**(-2.0_r2 * belt_alpha(k)) + s**(-2.0_r2 * belt_alpha_out(k)))

    case default
      prof = 0.0_r2

    end select

  end function radial_profile_fn


  pure function density_belt(X, Y, Zd, k) result(rho)       ! one belt, without normalisation

    real(r2), intent(in) :: X, Y, Zd                         ! disk-frame coordinates [m]
    integer,  intent(in) :: k
    real(r2)             :: rho
    real(r2)             :: r, varpi, s, cosf, a, h

    rho   = 0.0_r2
    varpi = sqrt(X*X + Y*Y)
    r     = sqrt(varpi*varpi + Zd*Zd)
    if (r <= 0.0_r2) return

    if (switch_vertical == 1) then
      if (abs(Zd) > sin_open * r) return
      s = r
    else
      if (varpi <= 0.0_r2) return
      if (abs(Zd) > tan_open * varpi) return
      s = varpi
    end if

    if (s < belt_rin_eff(k) * (1.0_r2 - 1.0e-12_r2) .or. s > belt_rout_eff(k) * (1.0_r2 + 1.0e-12_r2)) return

    if (belt_ecc(k) > 0.0_r2) then
      if (varpi > 0.0_r2) then
        cosf = (Y * belt_cos_om(k) + X * belt_sin_om(k)) / varpi     ! cos(psi - omega), psi = atan2(X, Y)
      else
        cosf = 0.0_r2
      end if
      a = s * (1.0_r2 + belt_ecc(k) * cosf) / (1.0_r2 - belt_ecc(k)**2)
    else
      a = s
    end if

    if (a < belt_Rmin(k) .or. a > belt_Rmax(k)) return

    if (belt_law(k) == 5) then
      rho = user_density_fn(a, X, Y, Zd, k, user_par)
    else
      rho = radial_profile_fn(a, k)
    end if
    rho = rho * (a / s)

    if (switch_vertical == 2) then
      h   = h_over_r * varpi
      rho = rho * exp(-0.5_r2 * (Zd / h)**2)
    end if

  end function density_belt


  pure function density_disk(X, Y, Zd) result(rho)

    real(r2), intent(in) :: X, Y, Zd                         ! disk-frame coordinates [m]
    real(r2)             :: rho
    integer              :: k

    if (n_belts == 1) then
      rho = belt_amp(1) * density_belt(X, Y, Zd, 1)
    else
      rho = 0.0_r2
      do k = 1, n_belts
        rho = rho + belt_amp(k) * density_belt(X, Y, Zd, k)
      end do
    end if

  end function density_disk


  pure function density_fn(xd, yd, z) result(rho)          ! density at a detector-frame point

    real(r2), intent(in) :: xd, yd, z
    real(r2)             :: rho
    real(r2)             :: X, Y, Zd

    call detector_to_disk(xd, yd, z, X, Y, Zd)
    rho = density_disk(X, Y, Zd)

  end function density_fn

  !_______________________________________________________________________________________|
  ! 3. Planck function and stellar intensity                                               |
  !_______________________________________________________________________________________|

  pure function Planck_fn(temperature, lambda) result(planck_result)

    real(r2), intent(in) :: temperature                      ! [K]
    real(r2), intent(in) :: lambda                           ! [m]
    real(r2)             :: planck_result                    ! B_lambda [W/m^2/m/sr]
    real(r2)             :: xexp

    planck_result = 0.0_r2
    if (temperature <= 0.0_r2 .or. lambda <= 0.0_r2) return

    xexp = (h * c) / (k_B * lambda * temperature)

    if (xexp > 700.0_r2) then                                ! exp() would overflow: B ~ 0
      planck_result = 0.0_r2
    else if (xexp < 1.0e-6_r2) then                          ! Rayleigh-Jeans limit
      planck_result = 2.0_r2 * c * k_B * temperature / lambda**4
    else
      planck_result = (2.0_r2 * h * c * c / lambda**5) / (exp(xexp) - 1.0_r2)
    end if

  end function Planck_fn


  pure function stellar_radiance(lambda) result(Istar)      ! surface intensity [W/m^2/m/sr]

    real(r2), intent(in) :: lambda                           ! [m]
    real(r2)             :: Istar
    integer              :: i
    real(r2)             :: t

    if (switch_BB_or_readingspec == 1) then

      Istar = Planck_fn(Tstar, lambda)

    else

      if (lambda < spec_lambda(1)) then
        Istar = 0.0_r2
      else if (lambda >= spec_lambda(n_spec)) then           ! Rayleigh-Jeans continuation
        Istar = spec_I(n_spec) * (spec_lambda(n_spec) / lambda)**4
      else
        i = locate(spec_lambda, lambda)
        if (spec_I(i) > 0.0_r2 .and. spec_I(i+1) > 0.0_r2) then
          t     = log(lambda / spec_lambda(i)) / log(spec_lambda(i+1) / spec_lambda(i))
          Istar = exp(log(spec_I(i)) + t * log(spec_I(i+1) / spec_I(i)))
        else
          t     = (lambda - spec_lambda(i)) / (spec_lambda(i+1) - spec_lambda(i))
          Istar = spec_I(i) + t * (spec_I(i+1) - spec_I(i))
        end if
      end if

    end if

  end function stellar_radiance

  !_______________________________________________________________________________________|
  ! 4. Henyey-Greenstein phase function, normalised to int p dOmega = 1                    |
  !_______________________________________________________________________________________|

  pure function hg_phase(g, cos_theta) result(p)

    real(r2), intent(in) :: g, cos_theta
    real(r2)             :: p

    p = (1.0_r2 - g*g) / (4.0_r2 * pi * (1.0_r2 + g*g - 2.0_r2 * g * cos_theta)**1.5_r2)

  end function hg_phase

  !_______________________________________________________________________________________|
  ! 5. Look-up helpers                                                                     |
  !_______________________________________________________________________________________|

  pure function dist_index(r) result(id)                    ! nearest point of the (logarithmic)
                                                            ! distance grid
    real(r2), intent(in) :: r
    integer              :: id

    if (r <= dist_lo) then
      id = 1
      return
    end if
    id = nint(log(r / dist_lo) / dist_step) + 1
    if (id < 1)      id = 1
    if (id > n_dist) id = n_dist

  end function dist_index


  pure function locate(xarr, x) result(i)                   ! xarr(i) <= x < xarr(i+1), monotonic

    real(r2), dimension(:), intent(in) :: xarr
    real(r2), intent(in)               :: x
    integer                            :: i, lo, hi, mid, n
    logical                            :: ascending

    n = size(xarr)
    ascending = (xarr(n) >= xarr(1))
    lo = 1
    hi = n

    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if ((x >= xarr(mid)) .eqv. ascending) then
        lo = mid
      else
        hi = mid
      end if
    end do

    i = lo
    if (i > n - 1) i = n - 1
    if (i < 1)     i = 1

  end function locate


  pure function interp_lin(xarr, yarr, x) result(y)         ! linear interpolation, clamped ends

    real(r2), dimension(:), intent(in) :: xarr, yarr
    real(r2), intent(in)               :: x
    real(r2)                           :: y
    integer                            :: i
    real(r2)                           :: t

    i = locate(xarr, x)
    t = (x - xarr(i)) / (xarr(i+1) - xarr(i))
    t = max(0.0_r2, min(1.0_r2, t))
    y = yarr(i) + t * (yarr(i+1) - yarr(i))

  end function interp_lin


  pure subroutine box_distances(cx, cy, cz, hh, dmin, dmax) ! min./max. distance of a cube to (0,0,0)

    real(r2), intent(in)  :: cx, cy, cz, hh
    real(r2), intent(out) :: dmin, dmax

    dmin = sqrt(max(abs(cx) - hh, 0.0_r2)**2 + max(abs(cy) - hh, 0.0_r2)**2 + &
                max(abs(cz) - hh, 0.0_r2)**2)
    dmax = sqrt((abs(cx) + hh)**2 + (abs(cy) + hh)**2 + (abs(cz) + hh)**2)

  end subroutine box_distances


  pure function cube_outside_wedge(cx, cy, cz, hh) result(outside)

    ! .true. if no point of the cube (centre cx,cy,cz in the detector frame, half size hh)
    ! can lie inside the wedge |Z| <= sin(theta) r. Conservative test.

    real(r2), intent(in) :: cx, cy, cz, hh
    logical              :: outside
    real(r2)             :: X, Y, Zd, rc, hdiag

    call detector_to_disk(cx, cy, cz, X, Y, Zd)
    rc      = sqrt(cx*cx + cy*cy + cz*cz)
    hdiag   = sqrt(3.0_r2) * hh
    outside = (abs(Zd) - hdiag) > sin_open * (rc + hdiag)

  end function cube_outside_wedge

  !_______________________________________________________________________________________|
  ! 6. Robust reader for numerical tables (lines that cannot be read are skipped,          |
  !    e.g. headers, comments or a leading "number of lines" entry)                        |
  !_______________________________________________________________________________________|

  subroutine read_table(fname, ncol, tab, nrow, min_rows)

    character(len=*), intent(in)                        :: fname
    integer, intent(in)                                 :: ncol
    real(r2), allocatable, dimension(:,:), intent(out)  :: tab        ! (ncol, nrow)
    integer, intent(out)                                :: nrow
    integer, intent(in), optional                       :: min_rows   ! default 2
    integer                                             :: nmin

    character(len=1024)                   :: line
    real(r2), dimension(ncol)             :: vals
    real(r2), allocatable, dimension(:,:) :: tmp
    integer                               :: unit, ios, cap
    logical                               :: ex, isdir

    nmin = 2
    if (present(min_rows)) nmin = min_rows

    isdir = .false.
    inquire(file=fname, exist=ex)
    if (ex) inquire(file=trim(fname)//"/.", exist=isdir)
    if (ex .and. isdir) then
      call dms_fail("file name missing (directory given): "//trim(fname))
      nrow = 0
      allocate(tab(ncol, 0))
      return
    end if
    if (.not. ex) then
      call dms_fail("file not found: "//trim(fname))
      nrow = 0
      allocate(tab(ncol, 0))
      return
    end if

    cap  = 1024
    allocate(tmp(ncol, cap))
    nrow = 0

    open(newunit=unit, file=fname, action="read", status="old", form="formatted")

    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '%') cycle
      read(line, *, iostat=ios) vals
      if (ios /= 0) cycle
      if (nrow == cap) then
        call grow(tmp, cap)
      end if
      nrow = nrow + 1
      tmp(:, nrow) = vals
    end do

    close(unit)

    if (nrow < nmin) then
      call dms_fail("too few valid data lines ("//integer_to_string(nrow)//") in "//trim(fname))
      deallocate(tmp)
      allocate(tab(ncol, 0))
      return
    end if

    allocate(tab(ncol, nrow))
    tab(:, :) = tmp(:, 1:nrow)
    deallocate(tmp)

  contains

    subroutine grow(arr, n)
      real(r2), allocatable, dimension(:,:), intent(inout) :: arr
      integer, intent(inout)                               :: n
      real(r2), allocatable, dimension(:,:)                :: t2
      allocate(t2(size(arr,1), 2*n))
      t2(:, 1:n) = arr(:, 1:n)
      call move_alloc(t2, arr)
      n = 2 * n
    end subroutine grow

  end subroutine read_table

  !_______________________________________________________________________________________|
  ! 7. Strings and file names                                                              |
  !_______________________________________________________________________________________|

  function integer_to_string(i) result(result_itos)

    character(:), allocatable :: result_itos
    integer, intent(in)       :: i
    character(len=24)         :: tmp1

    write(tmp1, '(i0)') i
    result_itos = trim(adjustl(tmp1))

  end function integer_to_string


  function real_to_string(x) result(result_rtos)

    character(:), allocatable :: result_rtos
    real(r2), intent(in)      :: x
    character(len=40)         :: tmp2

    if (abs(x) >= 1.0e7_r2 .or. (abs(x) < 1.0e-2_r2 .and. x /= 0.0_r2)) then
      if (abs(x) < 1.0e-99_r2 .or. abs(x) >= 1.0e100_r2) then
        write(tmp2, '(es11.3e3)') x
      else
        write(tmp2, '(es10.3)') x
      end if
    else
      write(tmp2, '(f0.2)') x
      if (tmp2(1:1) == '.')  tmp2 = '0'//tmp2(1:39)
      if (tmp2(1:2) == '-.') tmp2 = '-0'//tmp2(2:39)
    end if
    result_rtos = trim(adjustl(tmp2))

  end function real_to_string


  function real_to_scientific_string(x) result(result_rtoes)

    character(:), allocatable :: result_rtoes
    real(r2), intent(in)      :: x
    character(len=40)         :: tmp3

    if ((abs(x) < 1.0e-99_r2 .and. x /= 0.0_r2) .or. abs(x) >= 1.0e100_r2) then
      write(tmp3, '(es13.5e3)') x
    else
      write(tmp3, '(es12.5)') x
    end if
    result_rtoes = trim(adjustl(tmp3))

  end function real_to_scientific_string


  function file_suffix(il) result(suf)

    ! File name suffix following switch_file_extensions:
    !   0: "_#il"                         1: ".name_#il_Rin_.._alpha_.._q_.._amin_.._e_.._inc_.."
    !   2: ".name_#il"
    ! il <= 0 drops the wavelength tag (used for SEDs and tables).

    integer, intent(in)       :: il
    character(:), allocatable :: suf

    suf = ""
    if (switch_file_extensions == 1 .or. switch_file_extensions == 2) then
      suf = "."//trim(adjustl(file_extension_name))
    end if
    if (il > 0) suf = suf//"_#"//integer_to_string(il)
    if (switch_file_extensions == 1) then
      suf = suf//"_Rin_"//real_to_string(Rmin_au)//"_alpha_"//real_to_string(alpha) &
               //"_q_"//real_to_string(q_size)//"_amin_"//real_to_string(grain_radius_min) &
               //"_e_"//real_to_string(ecc)//"_inc_"//real_to_string(inc_deg)
    end if

  end function file_suffix

end module functions
