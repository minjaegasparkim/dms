subroutine dust_mixing

  !___________________________________________________________________________________________________|
  !                                                                                                   |
  ! Effective medium theory                                                                           |
  !___________________________________________________________________________________________________|
  ! The effective dielectric function of a mixture of a matrix (eps_m) with inclusions (eps_i,         |
  ! volume fraction f) is computed with one of four mixing rules. The result is written to            |
  ! ../ri-data/<emt_output_file> (columns: wavelength [micron], n, k) and can then be used as a       |
  ! dust component like any other optical-data file.                                                  |
  !                                                                                                   |
  !  Maxwell Garnett : eps = eps_m [1 + 3 f beta / (1 - f beta)],  beta = (eps_i-eps_m)/(eps_i+2eps_m)  |
  !  Bruggeman       : f (eps_i-eps)/(eps_i+2eps) + (1-f)(eps_m-eps)/(eps_m+2eps) = 0                  |
  !                    -> eps = [b + sqrt(b^2 + 8 eps_i eps_m)]/4,  b = (3f-1) eps_i + (2-3f) eps_m     |
  !  Lichtenecker    : eps = eps_i^f eps_m^(1-f)                                                      |
  !  Looyenga        : eps = [f eps_i^(1/3) + (1-f) eps_m^(1/3)]^3                                    |
  !                                                                                                   |
  ! References: Maxwell Garnett (1904, Phil. Trans. R. Soc. Lond. 203, 385); Bruggeman (1935,         |
  ! Ann. Phys. 24, 636); Looyenga (1965, Physica 31, 401); Maron & Maron (2005, MNRAS 357, 873)       |
  !                                                                                                   |
  ! Changes: the matrix and inclusion files may now have different wavelength grids (the inclusion    |
  ! data are interpolated onto the matrix grid within the common range); the complex arithmetic is    |
  ! done in double precision (before: default single-precision COMPLEX); the physical Bruggeman root  |
  ! (Im(eps) >= 0) is selected explicitly.                                                            |
  !---------------------------------------------------------------------------------------------------|

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2), allocatable, dimension(:,:) :: tab_m, tab_i
  integer     :: n_m, n_i, il, unit
  real(r2)    :: lam, n_in, k_in, nmix, kmix, lam_lo, lam_hi
  complex(r2) :: e_m, e_in, e_av, beta, b, sq, e1, e2

  if (switch_dust_mix /= 1) return

  if (verbose) print *, " "
  if (verbose) print *, " 0) effective medium approximation"
  if (verbose) print *, " "

  if (vol_frac <= 0.0_r2 .or. vol_frac >= 1.0_r2) then
    call dms_fail("the volume fraction of the inclusion must be in ]0, 1[")
    return
  end if

  call read_table(trim(ri_dir)//trim(emt_matrix_file),    3, tab_m, n_m)
  call read_table(trim(ri_dir)//trim(emt_inclusion_file), 3, tab_i, n_i)
  if (dms_ierr /= 0) return

  lam_lo = max(minval(tab_m(1,:)), minval(tab_i(1,:)))
  lam_hi = min(maxval(tab_m(1,:)), maxval(tab_i(1,:)))

  if (verbose) print *, "   matrix    : ", trim(emt_matrix_file),    "  (", real_to_string((1.0_r2-vol_frac)*100.0_r2), " %)"
  if (verbose) print *, "   inclusion : ", trim(emt_inclusion_file), "  (", real_to_string(vol_frac*100.0_r2), " %)"

  open(newunit=unit, file=trim(ri_dir)//trim(emt_output_file), action="write", status="replace", &
       form="formatted")
  write(unit, '(a)') "# EMT mixture written by DMS: wavelength [micron], n, k"

  do il = 1, n_m

    lam = tab_m(1, il)
    if (lam < lam_lo .or. lam > lam_hi) cycle

    n_in = interp_lin(log(tab_i(1,:)), tab_i(2,:), log(lam))
    k_in = interp_lin(log(tab_i(1,:)), tab_i(3,:), log(lam))

    e_m  = cmplx(tab_m(2,il)**2 - tab_m(3,il)**2, 2.0_r2 * tab_m(2,il) * tab_m(3,il), kind=r2)
    e_in = cmplx(n_in**2 - k_in**2, 2.0_r2 * n_in * k_in, kind=r2)

    select case (switch_dust_mix_method)

    case (1)                                                   ! Maxwell Garnett
      beta = (e_in - e_m) / (e_in + 2.0_r2 * e_m)
      e_av = e_m * (1.0_r2 + 3.0_r2 * vol_frac * beta / (1.0_r2 - vol_frac * beta))

    case (2)                                                   ! Bruggeman
      b  = (3.0_r2 * vol_frac - 1.0_r2) * e_in + (2.0_r2 - 3.0_r2 * vol_frac) * e_m
      sq = sqrt(b * b + 8.0_r2 * e_in * e_m)
      e1 = (b + sq) / 4.0_r2
      e2 = (b - sq) / 4.0_r2
      if (aimag(e1) >= 0.0_r2 .and. aimag(e2) >= 0.0_r2) then
        if (real(e1) >= real(e2)) then
          e_av = e1
        else
          e_av = e2
        end if
      else if (aimag(e1) >= aimag(e2)) then
        e_av = e1
      else
        e_av = e2
      end if

    case (3)                                                   ! Lichtenecker
      e_av = e_in**vol_frac * e_m**(1.0_r2 - vol_frac)

    case (4)                                                   ! Looyenga
      e_av = (vol_frac * e_in**(1.0_r2/3.0_r2) + (1.0_r2 - vol_frac) * e_m**(1.0_r2/3.0_r2))**3

    end select

    nmix = sqrt(max(0.0_r2, (abs(e_av) + real(e_av)) / 2.0_r2))
    kmix = sqrt(max(0.0_r2, (abs(e_av) - real(e_av)) / 2.0_r2))

    write(unit, '(3es17.8e3)') lam, nmix, kmix

  end do

  close(unit)

  if (verbose) print *, "   mixture written to ", trim(ri_dir)//trim(emt_output_file)
  if (verbose) print *, " "

end subroutine dust_mixing
