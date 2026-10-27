!____________________________________________________________________________________________________
!
! 10th module : output of images, radial profiles and SEDs (ASCII and/or FITS)
!____________________________________________________________________________________________________
!
! Images
!   ASCII : ../results/<name><suffix>.dat, one file per wavelength (suffix "_#il"):
!           a header (lines starting with '#'), then N rows x N columns;
!           row j = detector y index iy = j (yd increasing), column i = detector x index ix.
!           numpy: np.loadtxt(file)  ->  array[iy-1, ix-1]
!   FITS  : ../results/<name><suffix>.fits, one cube (NAXIS1 = x, NAXIS2 = y, NAXIS3 = wavelength)
!           with the wavelengths in the header (WAVE0001 ...), pixel scale in arcsec and au,
!           and the main model parameters.
!   Pixel (ix, iy) is centred at xd = (ix - (N+1)/2) * pixel size, yd = (iy - (N+1)/2) * pixel size;
!   the star is in the centre of pixel ((N+1)/2, (N+1)/2).
!
! SED
!   ASCII : ../results/SED<suffix>.dat  (columns described in the header)
!   FITS  : ../results/SED<suffix>.fits (2-d array: wavelength x column), optional
!____________________________________________________________________________________________________

module output_routines

  use datatype
  use constants
  use parameters
  use functions
  use fits_io

  implicit none

contains

  !-----------------------------------------------------------------------------------------------
  ! common FITS header cards describing the model
  !-----------------------------------------------------------------------------------------------
  subroutine model_cards(hdr)

    type(fits_header), intent(inout) :: hdr

    call fits_card_str (hdr, 'ORIGIN',   'DMS (Debris disk Modelling Software)', '')
    call fits_card_real(hdr, 'DIST_PC',  dist_star_pc,       'distance [pc]')
    call fits_card_real(hdr, 'TSTAR',    Tstar,              'stellar effective temperature [K]')
    call fits_card_real(hdr, 'RSTAR',    Rstar_Rsun,         'stellar radius [R_sun]')
    call fits_card_real(hdr, 'MSTAR',    Mstar_Msun,         'stellar mass [M_sun]')
    call fits_card_real(hdr, 'LSTAR',    Lstar / 3.828e26_r2, 'stellar luminosity [L_sun]')
    call fits_card_real(hdr, 'RMIN_AU',  Rmin_au,            'inner semi-major axis [au]')
    call fits_card_real(hdr, 'RMAX_AU',  Rmax_au,            'outer semi-major axis [au]')
    call fits_card_real(hdr, 'ECC',      ecc,                'eccentricity')
    call fits_card_real(hdr, 'OMEGA',    omega_peri_deg,     'argument of pericentre [deg]')
    call fits_card_real(hdr, 'INC',      inc_deg,            'inclination [deg]')
    call fits_card_real(hdr, 'PA',       pa_deg,             'position angle [deg]')
    call fits_card_real(hdr, 'OPENANG',  open_angle_deg,     'half opening angle [deg]')
    call fits_card_int (hdr, 'DENSLAW',  switch_density_fn,  'density law')
    call fits_card_real(hdr, 'ALPHA',    alpha,              'density exponent')
    call fits_card_real(hdr, 'MDUST',    disk_mass_Msun,     'dust mass [M_sun]')
    call fits_card_real(hdr, 'AMIN',     grain_radius_min,   'minimum grain radius [micron]')
    call fits_card_real(hdr, 'AMAX',     grain_radius_max,   'maximum grain radius [micron]')
    call fits_card_real(hdr, 'QSIZE',    q_size,             'size distribution exponent')
    call fits_card_int (hdr, 'PHASEFN',  switch_phase_function, '1 Mie, 2 HG, 3 two-component HG')
    call fits_card_int (hdr, 'SUBLIM',   switch_sublimation, 'sublimation 0/1')
    call fits_card_int (hdr, 'BLOWOUT',  switch_blowout,     'blow-out removal 0/1')
    if (switch_blowout == 1) then
      call fits_card_real(hdr, 'BETACMIN', beta_crit_min,    'smallest blow-out threshold beta')
      call fits_card_real(hdr, 'BETACMAX', beta_crit_max,    'largest blow-out threshold beta')
    end if
    call fits_card_real(hdr, 'FDUST',    f_dust,             'fractional luminosity L_d/L*')

  end subroutine model_cards

  !-----------------------------------------------------------------------------------------------
  ! one image quantity, all wavelengths: img(nz, N, N) with nz = n_wavelength_obs or 1
  !-----------------------------------------------------------------------------------------------
  subroutine write_image(name, img, unit_str, comment, per_wavelength)

    character(len=*), intent(in)              :: name, unit_str, comment
    real(r2), dimension(:,:,:), intent(in)    :: img          ! (nz, ix, iy)
    logical, intent(in)                       :: per_wavelength   ! .false.: maps (tau, mass)

    integer  :: nz, n, iz, ix, iy, u
    real(r2), allocatable, dimension(:,:,:) :: cube
    type(fits_header) :: hdr
    character(len=8)  :: key
    character(len=:), allocatable :: fname

    nz = size(img, 1)
    n  = size(img, 2)

    if (switch_output_format == 1 .or. switch_output_format == 3) then
      do iz = 1, nz
        if (.not. per_wavelength) then
          fname = trim(results_dir)//name//file_suffix(0)//".dat"
        else
          fname = trim(results_dir)//name//file_suffix(iz)//".dat"
        end if
        open(newunit=u, file=fname, action="write", status="replace", form="formatted")
        write(u, '(a)') "# DMS image: "//comment
        if (per_wavelength) write(u, '(a, es15.6e3, a)') "# wavelength = ", wavelength_obs(iz), " micron"
        write(u, '(a)') "# unit: "//unit_str
        write(u, '(a, i0, a, i0, a, es15.6e3, a, es15.6e3, a)') "# ", n, " x ", n, " pixels; pixel = ", &
              pixel_scale_arcsec, " arcsec = ", pixel_size_m / au, " au"
        write(u, '(a)') "# rows: y index 1..N (yd increasing), columns: x index 1..N; star at the centre"
        do iy = 1, n
          write(u, '(*(es14.5e3))') (img(iz, ix, iy), ix = 1, n)
        end do
        close(u)
      end do
    end if

    if (switch_output_format == 2 .or. switch_output_format == 3) then
      allocate(cube(n, n, nz))
      do iz = 1, nz
        cube(:, :, iz) = img(iz, :, :)
      end do
      call fits_header_init(hdr)
      call fits_card_str (hdr, 'BUNIT',  unit_str, '')
      call fits_card_str (hdr, 'CONTENT', name, comment)
      call fits_card_str (hdr, 'CTYPE1', 'OFFSET-X', 'detector x offset from the star')
      call fits_card_str (hdr, 'CUNIT1', 'arcsec', '')
      call fits_card_real(hdr, 'CRPIX1', real((n + 1) / 2, kind=r2), 'star pixel')
      call fits_card_real(hdr, 'CRVAL1', 0.0_r2, '')
      call fits_card_real(hdr, 'CDELT1', pixel_scale_arcsec, 'pixel scale [arcsec]')
      call fits_card_str (hdr, 'CTYPE2', 'OFFSET-Y', 'detector y offset from the star')
      call fits_card_str (hdr, 'CUNIT2', 'arcsec', '')
      call fits_card_real(hdr, 'CRPIX2', real((n + 1) / 2, kind=r2), 'star pixel')
      call fits_card_real(hdr, 'CRVAL2', 0.0_r2, '')
      call fits_card_real(hdr, 'CDELT2', pixel_scale_arcsec, 'pixel scale [arcsec]')
      call fits_card_real(hdr, 'PIXAU',  pixel_size_m / au, 'pixel size [au]')
      if (per_wavelength) then
        call fits_card_int(hdr, 'NWAVE', nz, 'number of wavelengths (axis 3)')
        do iz = 1, min(nz, 9999)
          write(key, '(a, i4.4)') 'WAVE', iz
          call fits_card_real(hdr, key, wavelength_obs(iz), 'wavelength [micron]')
        end do
      end if
      call model_cards(hdr)
      fname = trim(results_dir)//name//file_suffix(0)//".fits"
      call fits_write_cube(fname, cube, hdr)
      deallocate(cube)
    end if

  end subroutine write_image

  !-----------------------------------------------------------------------------------------------
  ! radial brightness profiles along the projected major axis (disk line of nodes),
  ! bilinear interpolation of the images; direction (-sin PA, cos PA) in the detector frame
  !-----------------------------------------------------------------------------------------------
  subroutine write_profiles(n, ufac, ustr)

    integer, intent(in)          :: n
    real(r2), intent(in)         :: ufac          ! unit conversion from Jy/pixel
    character(len=*), intent(in) :: ustr
    integer  :: il, k, u, nq, ic0
    real(r2) :: s, xd, yd
    real(r2), allocatable :: vals(:)
    character(len=:), allocatable :: fname, head

    ic0 = (n + 1) / 2
    head = "# 1 offset [au]  2 offset [arcsec]"
    nq = 0
    if (do_reem) then
      nq = nq + 1; head = head//"  "//integer_to_string(nq+2)//" thermal"
    end if
    if (do_sca) then
      head = head//"  "//integer_to_string(nq+3)//" scattered  "//integer_to_string(nq+4)//" Q  "// &
             integer_to_string(nq+5)//" U  "//integer_to_string(nq+6)//" PI  "//integer_to_string(nq+7)//" Q_phi"
      nq = nq + 5
    end if
    nq = nq + 1; head = head//"  "//integer_to_string(nq+2)//" total   ["//ustr//"]"
    allocate(vals(nq))

    do il = 1, n_wavelength_obs
      fname = trim(results_dir)//"profile_major_axis"//file_suffix(il)//".dat"
      open(newunit=u, file=fname, action="write", status="replace", form="formatted")
      write(u, '(a, es15.6e3, a)') "# radial profile along the projected major axis, wavelength = ", &
            wavelength_obs(il), " micron"
      write(u, '(a)') "# positive offsets towards (xd, yd) = (-sin PA, cos PA); bilinear interpolation"
      write(u, '(a)') head
      do k = -(ic0 - 1), ic0 - 1
        s  = real(k, kind=r2) * pixel_size_m
        xd = -s * sin_pa
        yd =  s * cos_pa
        nq = 0
        if (do_reem) then
          nq = nq + 1; vals(nq) = bilinear(image_reem, il, xd, yd, n)
        end if
        if (do_sca) then
          nq = nq + 1; vals(nq) = bilinear(image_sca,  il, xd, yd, n)
          nq = nq + 1; vals(nq) = bilinear(image_Q,    il, xd, yd, n)
          nq = nq + 1; vals(nq) = bilinear(image_U,    il, xd, yd, n)
          nq = nq + 1; vals(nq) = bilinear(image_pol,  il, xd, yd, n)
          nq = nq + 1; vals(nq) = bilinear(image_Qphi, il, xd, yd, n)
        end if
        nq = nq + 1; vals(nq) = bilinear(image_total, il, xd, yd, n)
        write(u, '(*(es15.6e3))') s / au, s / au / dist_star_pc, vals * ufac
      end do
      close(u)
    end do

    deallocate(vals)

  end subroutine write_profiles


  pure function bilinear(img, il, xd, yd, n) result(v)

    real(r2), dimension(:,:,:), intent(in) :: img
    integer, intent(in)                    :: il, n
    real(r2), intent(in)                   :: xd, yd
    real(r2) :: v, fx, fy, tx, ty
    integer  :: i0, j0

    fx = xd / pixel_size_m + real((n + 1) / 2, kind=r2)      ! fractional pixel index
    fy = yd / pixel_size_m + real((n + 1) / 2, kind=r2)
    i0 = max(1, min(n - 1, floor(fx)))
    j0 = max(1, min(n - 1, floor(fy)))
    tx = max(0.0_r2, min(1.0_r2, fx - real(i0, kind=r2)))
    ty = max(0.0_r2, min(1.0_r2, fy - real(j0, kind=r2)))
    v  = (1.0_r2 - tx) * (1.0_r2 - ty) * img(il, i0,   j0  ) + tx * (1.0_r2 - ty) * img(il, i0+1, j0  ) &
       + (1.0_r2 - tx) *          ty  * img(il, i0,   j0+1) + tx *          ty  * img(il, i0+1, j0+1)

  end function bilinear

  !-----------------------------------------------------------------------------------------------
  ! SED
  !-----------------------------------------------------------------------------------------------
  subroutine write_sed(origin)

    character(len=*), intent(in) :: origin
    integer :: il, u
    real(r2), allocatable, dimension(:,:,:) :: tab
    type(fits_header) :: hdr
    character(len=:), allocatable :: fname

    fname = trim(results_dir)//"SED"//file_suffix(0)//".dat"
    open(newunit=u, file=fname, action="write", status="replace", form="formatted")
    write(u, '(a)') "# DMS spectral energy distribution ("//origin//")"
    write(u, '(a, es15.6e3, a)') "# distance = ", dist_star_pc, " pc;  all fluxes in Jy"
    write(u, '(a)') "# 1 wavelength [micron]  2 thermal emission  3 scattered light  4 Stokes Q" // &
                    "  5 Stokes U  6 polarised flux sqrt(Q^2+U^2)  7 dust (2+3)  8 star  9 total (7+8)"
    do il = 1, n_wavelength_obs
      write(u, '(*(es15.6e3))') wavelength_obs(il), sed_reem(il), sed_sca(il), sed_Q(il), sed_U(il), &
                              sed_pol(il), sed_reem(il) + sed_sca(il), stellar_flux_Jy(il), &
                              sed_total(il)
    end do
    close(u)

    if (switch_output_format == 2 .or. switch_output_format == 3) then
      allocate(tab(n_wavelength_obs, 9, 1))
      tab(:, 1, 1) = wavelength_obs(:)
      tab(:, 2, 1) = sed_reem(:)
      tab(:, 3, 1) = sed_sca(:)
      tab(:, 4, 1) = sed_Q(:)
      tab(:, 5, 1) = sed_U(:)
      tab(:, 6, 1) = sed_pol(:)
      tab(:, 7, 1) = sed_reem(:) + sed_sca(:)
      tab(:, 8, 1) = stellar_flux_Jy(:)
      tab(:, 9, 1) = sed_total(:)
      call fits_header_init(hdr)
      call fits_card_str(hdr, 'BUNIT',   'Jy', 'columns 2-9; column 1 in micron')
      call fits_card_str(hdr, 'CONTENT', 'SED', origin)
      call fits_card_str(hdr, 'COL1', 'wavelength', 'micron')
      call fits_card_str(hdr, 'COL2', 'thermal emission', 'Jy')
      call fits_card_str(hdr, 'COL3', 'scattered light', 'Jy')
      call fits_card_str(hdr, 'COL4', 'Stokes Q', 'Jy')
      call fits_card_str(hdr, 'COL5', 'Stokes U', 'Jy')
      call fits_card_str(hdr, 'COL6', 'polarised flux', 'Jy')
      call fits_card_str(hdr, 'COL7', 'dust total', 'Jy')
      call fits_card_str(hdr, 'COL8', 'star', 'Jy')
      call fits_card_str(hdr, 'COL9', 'total', 'Jy')
      call model_cards(hdr)
      call fits_write_cube(trim(results_dir)//"SED"//file_suffix(0)//".fits", tab, hdr)
      deallocate(tab)
    end if

  end subroutine write_sed

  !-----------------------------------------------------------------------------------------------
  ! scattered and polarised flux as a function of the scattering angle (whole disk)
  !-----------------------------------------------------------------------------------------------
  subroutine write_phase()

    integer :: u, ib, il
    character(len=:), allocatable :: fname

    fname = trim(results_dir)//"phase_function"//file_suffix(0)//".dat"
    open(newunit=u, file=fname, action="write", status="replace", form="formatted")
    write(u, '(a)') "# scattered stellar light of the whole disk as a function of the scattering angle theta"
    write(u, '(a)') "# (angle between star->grain and grain->observer; 1-degree bins), flux per bin [Jy]"
    write(u, '(a)') "# polarised flux > 0: polarisation perpendicular to the scattering plane"
    write(u, '(a)') "# sum over all bins = scattered flux of the SED; empty bins: angles not seen at this inclination"
    write(u, '(a)', advance='no') "# 1 theta [deg]"
    do il = 1, n_wavelength_obs
      write(u, '(a, i0, a, es11.4e2, a, i0, a)', advance='no') "  ", 2*il, " F_sca(", wavelength_obs(il), &
            " um)  ", 2*il+1, " F_pol"
    end do
    write(u, *)
    do ib = 1, n_phase_bins
      write(u, '(f8.2)', advance='no') real(ib, kind=r2) - 0.5_r2
      do il = 1, n_wavelength_obs
        write(u, '(2es15.6e3)', advance='no') phase_sca(il, ib), phase_pol(il, ib)
      end do
      write(u, *)
    end do
    close(u)

  end subroutine write_phase

  !-----------------------------------------------------------------------------------------------
  ! 3-d cubes (x, y, z) of the flux emitted towards the observer by every volume element, one FITS
  ! file per wavelength and quantity, plus the scattering angle of every volume element
  !-----------------------------------------------------------------------------------------------
  subroutine write_cubes()

    integer :: il, n
    type(fits_header) :: hdr
    real(r2), allocatable :: tmp(:,:,:)

    n = number_of_pixels
    allocate(tmp(n, n, n))
    do il = 1, n_wavelength_obs
      call fits_header_init(hdr)
      call fits_card_str (hdr, 'BUNIT', 'Jy', 'flux per volume element towards the observer')
      call fits_card_str (hdr, 'CTYPE3', 'Z', 'line of sight (+z towards the observer)')
      call fits_card_real(hdr, 'PIXAU', pixel_size_m / au, 'volume element size [au]')
      call fits_card_real(hdr, 'CDELT1', pixel_scale_arcsec, 'pixel scale [arcsec]')
      call fits_card_real(hdr, 'WAVE', wavelength_obs(il), 'wavelength [micron]')
      call model_cards(hdr)
      if (do_reem) then
        tmp = cube_reem(il, :, :, :)
        call fits_write_cube(trim(results_dir)//"cube_reem"//file_suffix(il)//".fits", tmp, hdr)
      end if
      if (do_sca) then
        tmp = cube_sca(il, :, :, :)
        call fits_write_cube(trim(results_dir)//"cube_sca"//file_suffix(il)//".fits", tmp, hdr)
      end if
    end do
    call fits_header_init(hdr)
    call fits_card_str (hdr, 'BUNIT', 'deg', 'scattering angle at the centre of each element')
    call fits_card_real(hdr, 'PIXAU', pixel_size_m / au, 'volume element size [au]')
    call fits_write_cube(trim(results_dir)//"cube_scattering_angle"//file_suffix(0)//".fits", cube_angle, hdr)
    deallocate(tmp)

  end subroutine write_cubes

end module output_routines
