!====================================================================================================
!
!   DMS : Debris disks around Main-sequence Stars
!
!   Thermal re-emission, scattered light and polarisation of optically thin debris disks:
!   images (ASCII / FITS), SEDs, fractional luminosity and geometrical optical depth, including
!   sublimation, radiation-pressure blow-out, temperature-dependent optical constants,
!   (two-component) Henyey-Greenstein phase functions, several (eccentric) belts, user-defined
!   density distributions, model-atmosphere grids and simulated observations (PSF, mask).
!
!   Usage :  dms model.nml          input file (Fortran namelist, see examples/)
!            dms                    interactive questions (answers can be piped: dms < answers.txt)
!   Threads: export OMP_NUM_THREADS=8
!
!   M. Kim, S. Wolf, T. Loehne, F. Kirchschlager, A. V. Krivov (2018), A&A 618, A38
!
!====================================================================================================

program dms

  use parameters

  implicit none

  character(len=1024) :: arg
  integer :: nargs

  library_mode = .false.
  nargs = command_argument_count()

  if (nargs >= 1) then
    call get_command_argument(1, arg)
    if (trim(arg) == "-h" .or. trim(arg) == "--help") then
      print '(a)', "usage: dms [input.nml]   (without an argument: interactive questions)"
      stop
    end if
    call read_input_file(trim(arg))
  else
    call input_questions
  end if

  if (dms_ierr /= 0) error stop 1

  call dms_run

  if (dms_ierr /= 0) error stop 1

end program dms
