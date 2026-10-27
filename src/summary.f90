subroutine summary(run_time)

  !---------------------------------------------------------------------------------------------------
  ! Summary of the model and of the results, printed to the screen and written to
  ! ../results/summary<suffix>.txt
  !---------------------------------------------------------------------------------------------------

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  real(r2), intent(in) :: run_time               ! [s]

  integer :: u, ic, il, k
  character(len=:), allocatable :: fname, ext

  if (.not. (verbose .or. write_output)) return
  u = -1
  if (write_output) then
    call make_dir(results_dir)
    fname = trim(results_dir)//"summary"//file_suffix(0)//".txt"
    open(newunit=u, file=fname, action="write", status="replace", form="formatted")
  end if

  call say("______________________________________________________________________________")
  call say(" ")
  call say(" Summary")
  call say("______________________________________________________________________________")
  call say(" ")
  call say(" Star")
  call say("   T_eff = "//real_to_string(Tstar)//" K,  R* = "//real_to_string(Rstar_Rsun)// &
           " R_sun,  M* = "//real_to_string(Mstar_Msun)//" M_sun,  L* = "// &
           real_to_scientific_string(Lstar / 3.828e26_r2)//" L_sun,  d = "//real_to_string(dist_star_pc)//" pc")
  call say(" Disk")
  do k = 1, n_belts
    if (n_belts > 1) call say("   belt "//integer_to_string(k)//" (mass fraction "// &
                              real_to_string(100.0_r2 * belt_mass_frac(k))//" %)")
    call say("   density law "//trim(law_name(belt_law(k)))//", semi-major axes "// &
             real_to_string(belt_Rmin(k)/au)//" - "//real_to_string(belt_Rmax(k)/au)//" au,  e = "// &
             real_to_string(belt_ecc(k))//",  omega = "//real_to_string(belt_omega(k)/deg2rad)//" deg")
  end do
  call say("   dust from "//real_to_string(r_in_eff/au)//" au (pericentre) to "// &
           real_to_string(r_out_eff/au)//" au (apocentre)")
  if (switch_vertical == 1) then
    call say("   vertical: wedge, half opening angle "//real_to_string(open_angle_deg)//" deg")
  else
    call say("   vertical: Gaussian, h/r = "//real_to_string(h_over_r)//", cut at "// &
             real_to_string(gauss_cut)//" h")
  end if
  call say("   inclination "//real_to_string(inc_deg)//" deg,  position angle "//real_to_string(pa_deg)//" deg")
  call say("   dust mass "//real_to_scientific_string(disk_mass_Msun)//" M_sun = "// &
           real_to_scientific_string(disk_mass / 7.342e22_r2)//" M_moon")
  call say(" Dust")
  do ic = 1, n_grain_components
    if (switch_Tdep_optics == 1 .and. ic == Tdep_comp) then
      call say("   component "//integer_to_string(ic)//": temperature-dependent optical data ("// &
               integer_to_string(n_Tdep)//" sets)")
    end if
    call say("   component "//integer_to_string(ic)//": "//trim(dustfile_name(ic))//", rho = "// &
             real_to_string(particle_bulk_density(ic))//" g/cm^3, abundance "// &
             real_to_string(100.0_r2 * abun(ic))//" %")
    if (switch_sublimation == 1 .and. T_subl(ic) > 0.0_r2) then
      if (subl_to_comp(ic) == 0) then
        call say("      sublimation at "//real_to_string(T_subl(ic))//" K -> vanishes")
      else
        call say("      sublimation at "//real_to_string(T_subl(ic))//" K -> becomes component "// &
                 integer_to_string(subl_to_comp(ic)))
      end if
    end if
  end do
  if (switch_blowout == 1) then
    if (beta_crit_max - beta_crit_min > 1.0e-12_r2) then
      call say("   blow-out: grains with beta >= "//real_to_string(beta_crit_min)//" - "// &
               real_to_string(beta_crit_max)//" removed (threshold depends on the belt)")
    else
      call say("   blow-out: grains with beta >= "//real_to_string(beta_crit_min)//" removed")
    end if
  else
    call say("   blow-out: off")
  end if
  call say("   grain radii "//real_to_string(grain_radius_min)//" - "//real_to_string(grain_radius_max)// &
           " micron, q = "//real_to_string(q_size)//", "//integer_to_string(n_grain_radius)//" bins")
  select case (switch_phase_function)
  case (1)
    call say("   phase function: Mie theory")
  case (2)
    call say("   phase function: Henyey-Greenstein, g = "//real_to_string(hg_g1)// &
             ", p_max = "//real_to_string(pol_max))
  case (3)
    call say("   phase function: two-component HG, g1 = "//real_to_string(hg_g1)//", g2 = "// &
             real_to_string(hg_g2)//", w = "//real_to_string(hg_w)//", p_max = "//real_to_string(pol_max))
  end select

  call say(" ")
  call say(" Results")
  call say("   fractional luminosity           f_d = L_d/L* = "//real_to_scientific_string(f_dust))
  call say("   total geometrical cross section sum N pi a^2 = "// &
           real_to_scientific_string(total_cross_section / au**2)//" au^2")
  call say("   mean geometrical optical depth  (cross section / disk area) = "// &
           real_to_scientific_string(tau_mean))
  call say("   peak vertical optical depth (azimuthal mean)            = "// &
           real_to_scientific_string(tau_vertical_peak))
  call say("   largest radial optical depth in the mid-plane           = "// &
           real_to_scientific_string(tau_radial_max))
  if (tau_radial_max > 0.1_r2) then
    call say("   WARNING: radial optical depth > 0.1 - the optically thin assumption of DMS does not hold")
  end if
  if (switch_image_or_SED /= 2) then
    call say("   image: "//integer_to_string(number_of_pixels)//" x "//integer_to_string(number_of_pixels)// &
             " pixels, "//real_to_scientific_string(pixel_scale_arcsec)//" arcsec/pixel = "// &
             real_to_string(pixel_size_m / au)//" au/pixel")
    call say("   max. temperature difference in the finest sub-volume: "//real_to_string(maxTempDiff)//" K")
  end if
  call say(" ")
  call say("   wavelength [micron]   dust thermal [Jy]   scattered [Jy]   polarised [Jy]   star [Jy]")
  do il = 1, n_wavelength_obs
    block
      character(len=120) :: line
      write(line, '(3x, es15.6e3, 4es17.6e3)') wavelength_obs(il), sed_reem(il), sed_sca(il), sed_pol(il), &
                                           stellar_flux_Jy(il)
      call say(trim(line))
    end block
  end do

  call say(" ")
  if (.not. write_output) then
    call say(" (no files written)")
  else
  call say(" Output files in "//trim(results_dir))
  ext = file_suffix(0)
  if (switch_image_or_SED /= 2) then
    if (switch_output_format /= 2) call say("   image_*"//ext//"_#<wavelength index>.dat  (ASCII images)")
    if (switch_output_format /= 1) call say("   image_*"//ext//".fits  (FITS cubes)")
    call say("   image_tau"//ext//"  (line-of-sight geometrical optical depth)")
    if (have_obs) call say("   image_*_obs"//ext//"  (simulated observations: PSF / mask)")
    call say("   profile_major_axis"//ext//"_#<wavelength index>.dat")
  end if
  if (switch_image_or_SED /= 1) then
    if (switch_output_format == 1) then
      call say("   SED"//ext//".dat")
    else
      call say("   SED"//ext//".dat and .fits")
    end if
  end if
  call say("   optical_depth"//ext//".dat")
  if (switch_temperature_distribution == 1) call say("   Temperature_Distribution"//ext//".dat")
  if (switch_grainsize_distribution == 1)   call say("   grain_size_distribution"//ext//".dat")
  if (switch_miex_save == 1 .and. .not. optics_from_cache) &
    call say("   optical properties: "//trim(miex_dir)//trim(file_saved_miex_name)//".miex")
  end if
  call say(" ")
  call say(" Run time: "//real_to_string(run_time)//" s")
  call say("______________________________________________________________________________")

  if (u /= -1) close(u)

contains

  function law_name(l) result(nm)
    integer, intent(in) :: l
    character(len=24)   :: nm
    select case (l)
    case (1);  nm = "power law"
    case (2);  nm = "Ertel peak"
    case (3);  nm = "Gaussian"
    case (4);  nm = "broken power law"
    case (5);  nm = "user-defined"
    case default;  nm = "?"
    end select
  end function law_name

  subroutine say(line)
    character(len=*), intent(in) :: line
    if (verbose) print '(a)', line
    if (u /= -1) write(u, '(a)') line
  end subroutine say

end subroutine summary
