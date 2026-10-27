subroutine mie_cals

  !---------------------------------------------------------------------------------------------------
  ! Mie calculation on the wavelength grid of the optical data ("ri grid"):
  !   Qabs_ri(ilr, ig, is)  for the temperature calculation
  !   Qpr_ri (ilr, ig, is)  for the radiation-pressure (blow-out) parameter beta
  ! Only efficiencies are needed here, so the scattering amplitudes are not computed
  ! (doSA = .false.), which makes this step much faster than before. The amplitudes are computed
  ! at the observing wavelengths only (Qabs_Qsca_S11_S12_cals).
  !
  ! Optionally the results of a previous run are loaded from ../miexres/<name>.miex. The file is
  ! only used if grain sizes, wavelengths, scattering angles and optical constants are identical;
  ! otherwise everything is recomputed (with a warning).
  !
  ! The loop over grain sizes, species and wavelengths runs in parallel (OpenMP).
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions
  use mie_routines

  implicit none

  integer :: ntot, it, ilr, ig, is, ier
  real(r2) :: sizepar, qext, qsca, qabs, qbk, qpr, alb, gsc
  complex(r2) :: ri
  complex(r2), dimension(3) :: sdum1, sdum2

  if (verbose) print *, " "
  if (verbose) print *, " 3) Mie calculation on the optical-data grid"
  if (verbose) print *, " "

  allocate(Qabs_ri(n_wavelength_int, n_grain_radius, n_species), &
           Qpr_ri (n_wavelength_int, n_grain_radius, n_species))

  if (switch_precal == 1) then
    call load_optics_cache()
    if (optics_from_cache) then
      if (verbose) print *, "    optical properties loaded from ", trim(miex_dir)//trim(file_saved_miex_name)//".miex"
      return
    end if
  end if

  ntot = n_wavelength_int * n_grain_radius * n_species

  !$omp parallel do schedule(dynamic, 16) default(shared) &
  !$omp private(it, ilr, ig, is, sizepar, ri, qext, qsca, qabs, qbk, qpr, alb, gsc, ier, sdum1, sdum2)
  do it = 1, ntot

    ilr = mod(it - 1, n_wavelength_int) + 1
    ig  = mod((it - 1) / n_wavelength_int, n_grain_radius) + 1
    is  = (it - 1) / (n_wavelength_int * n_grain_radius) + 1

    sizepar = 2.0_r2 * pi * rad(ig) / wavelength_int(ilr)
    ri      = cmplx(n_ri(ilr, is), k_ri(ilr, is), kind=r2)

    if (switch_grain_shape == 1) then
      call dhs_efficiencies(ri, sizepar, dhs_fmax, dhs_nf, qext, qsca, qabs, qpr, gsc, ier)
    else
      call shexqnn2(ri, sizepar, qext, qsca, qabs, qbk, qpr, alb, gsc, ier, sdum1, sdum2, .false., 2)
    end if

    if (ier /= 0) then
      !$omp critical (mie_error)
      call dms_fail("Mie calculation failed: radius "//real_to_string(rad(ig))//" micron, wavelength "// &
                    real_to_string(wavelength_int(ilr))//" micron")
      !$omp end critical (mie_error)
      cycle
    end if

    ! negative/NaN absorption efficiencies (numerical noise for k ~ 0) are set to zero
    if (.not. (qabs > 0.0_r2)) qabs = 0.0_r2

    Qabs_ri(ilr, ig, is) = qabs
    Qpr_ri (ilr, ig, is) = qpr

  end do
  !$omp end parallel do
  if (dms_ierr /= 0) return

  if (verbose) print *, "    done: ", ntot, " Mie calculations"

end subroutine mie_cals


subroutine load_optics_cache()

  use datatype
  use parameters

  implicit none

  character(len=8)  :: magic
  integer           :: unit, ios, nlr, ng, ns, nl, na, has_s, shp, nfd
  real(r2)          :: fmx
  real(r2), allocatable, dimension(:) :: wl_r, rad_r, wl_o
  real(r2), allocatable, dimension(:,:) :: nr, kr
  logical :: ex, ok
  character(len=:), allocatable :: fname

  fname = trim(miex_dir)//trim(file_saved_miex_name)//".miex"
  optics_from_cache = .false.

  inquire(file=fname, exist=ex)
  if (.not. ex) then
    if (verbose) print *, "    warning: ", fname, " not found - computing the optical properties"
    return
  end if

  open(newunit=unit, file=fname, access='stream', form='unformatted', status='old', &
       action='read', iostat=ios)
  if (ios /= 0) return

  read(unit, iostat=ios) magic, nlr, ng, ns, nl, na, has_s, shp, nfd, fmx
  ok = (ios == 0 .and. magic == 'DMSMIE05' .and. nlr == n_wavelength_int .and. &
        shp == switch_grain_shape .and. &
        ng == n_grain_radius .and. ns == n_species .and. nl == n_wavelength_obs .and. &
        na == n_angle2)
  if (ok .and. shp == 1) ok = (nfd == dhs_nf .and. abs(fmx - dhs_fmax) <= 1.0e-12_r2)
  if (ok .and. do_sca) ok = (has_s == 1)

  if (ok) then
    allocate(wl_r(nlr), rad_r(ng), wl_o(nl), nr(nlr, ns), kr(nlr, ns))
    read(unit, iostat=ios) wl_r, rad_r, wl_o, nr, kr
    ok = (ios == 0)
    if (ok) ok = all(abs(wl_r - wavelength_int) <= 1.0e-10_r2 * wavelength_int) .and. &
                 all(abs(rad_r - rad) <= 1.0e-10_r2 * rad) .and. &
                 all(abs(wl_o - wavelength_obs) <= 1.0e-10_r2 * wavelength_obs) .and. &
                 all(abs(nr - n_ri) <= 1.0e-10_r2 * abs(n_ri) + 1.0e-14_r2) .and. &
                 all(abs(kr - k_ri) <= 1.0e-10_r2 * abs(k_ri) + 1.0e-14_r2)
  end if

  if (.not. ok) then
    if (verbose) print *, "    warning: ", fname, " does not match the current parameters - recomputing"
    close(unit)
    return
  end if

  allocate(Cabs(nl, ng, ns), Csca(nl, ng, ns))
  read(unit, iostat=ios) Qabs_ri, Qpr_ri, Cabs, Csca
  if (has_s == 1) then
    allocate(S11_norm(nl, na, ng, ns), S12_norm(nl, na, ng, ns))
    read(unit, iostat=ios) S11_norm, S12_norm
  end if
  close(unit)

  if (ios /= 0) then
    if (verbose) print *, "    warning: ", fname, " is incomplete - recomputing"
    deallocate(Cabs, Csca)
    if (allocated(S11_norm)) deallocate(S11_norm, S12_norm)
    return
  end if

  optics_from_cache = .true.

end subroutine load_optics_cache


subroutine save_optics_cache()

  use datatype
  use parameters
  use functions, only: make_dir
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char

  implicit none

  integer :: unit, has_s, ios
  character(len=:), allocatable :: fname, ftmp
  character(len=24) :: pidstr

  interface
    function c_rename(oldname, newname) bind(c, name='rename') result(r)
      import :: c_char, c_int
      character(kind=c_char), dimension(*), intent(in) :: oldname, newname
      integer(c_int) :: r
    end function c_rename
    function c_getpid() bind(c, name='getpid') result(p)
      import :: c_int
      integer(c_int) :: p
    end function c_getpid
  end interface

  call make_dir(miex_dir)
  fname = trim(miex_dir)//trim(file_saved_miex_name)//".miex"

  ! Written to a temporary file that is then renamed: several DMS processes running in parallel
  ! (batch runs) never read a partly written cache, and two processes saving the same name do
  ! not mix their data.
  write(pidstr, '(i0)') c_getpid()
  ftmp = fname//".tmp"//trim(pidstr)

  has_s = 0
  if (allocated(S11_norm)) has_s = 1

  open(newunit=unit, file=ftmp, access='stream', form='unformatted', status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    if (verbose) print *, "    warning: cannot write ", ftmp, " - optical properties not saved"
    return
  end if
  write(unit) 'DMSMIE05', n_wavelength_int, n_grain_radius, n_species, n_wavelength_obs, &
              n_angle2, has_s, switch_grain_shape, dhs_nf, dhs_fmax
  write(unit) wavelength_int, rad, wavelength_obs, n_ri, k_ri
  write(unit) Qabs_ri, Qpr_ri, Cabs, Csca
  if (has_s == 1) write(unit) S11_norm, S12_norm
  close(unit)

  if (c_rename(ftmp//c_null_char, fname//c_null_char) /= 0) then
    if (verbose) print *, "    warning: cannot rename ", ftmp, " to ", fname
    open(newunit=unit, file=ftmp, status='old', iostat=ios)
    if (ios == 0) close(unit, status='delete')
    return
  end if

  if (verbose) print *, "    optical properties saved to ", fname

end subroutine save_optics_cache
