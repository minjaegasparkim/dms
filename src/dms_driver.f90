!____________________________________________________________________________________________________
!
! Driver: the sequence of a DMS run after the input has been read (used by the program dms and
! by the Python interface). After every step the error flag is checked; in library mode an
! error ends the run and is reported to the caller.
!____________________________________________________________________________________________________

subroutine dms_run

  use datatype
  use parameters
  use functions

  implicit none

  integer(i8) :: t_start, t_end, rate

  call system_clock(t_start, rate)
  run_complete = .false.

  call finalize_input;                                  if (dms_ierr /= 0) return
  if (switch_dust_mix == 1) call dust_mixing;           if (dms_ierr /= 0) return
  call pre_cals;                                        if (dms_ierr /= 0) return
  call scattering_angles_cals;                          if (dms_ierr /= 0) return
  call reading_spectrum;                                if (dms_ierr /= 0) return
  call reading_temperature_dependent_data;              if (dms_ierr /= 0) return
  call mie_cals;                                        if (dms_ierr /= 0) return
  call Qabs_Qsca_S11_S12_cals;                          if (dms_ierr /= 0) return
  call temperature_distribution_cals;                   if (dms_ierr /= 0) return
  call particle_abundances;                             if (dms_ierr /= 0) return
  call belt_normalisation;                              if (dms_ierr /= 0) return

  if (switch_image_or_SED == 1 .or. switch_image_or_SED == 3) then
    call max_temp_diff_cals;                            if (dms_ierr /= 0) return
    call computing_images;                              if (dms_ierr /= 0) return
    call observation_products;                          if (dms_ierr /= 0) return
  else
    call computing_SED;                                 if (dms_ierr /= 0) return
  end if

  call optical_depth_cals;                              if (dms_ierr /= 0) return

  call system_clock(t_end)
  call summary(real(t_end - t_start, kind=r2) / real(rate, kind=r2))

  run_complete = .true.

end subroutine dms_run
