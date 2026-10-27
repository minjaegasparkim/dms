!____________________________________________________________________________________________________
!
! 1. 1st module : datatype
!____________________________________________________________________________________________________
!================================================================================================
! Kind parameters used in the whole code.
!
!   r1 : default single precision real  (selected_real_kind(6))
!   r2 : double precision real          (selected_real_kind(15, 307))
!   i8 : 64-bit integer (used for array sizes and byte counters, e.g. in the FITS writer)
!
! (The previous version used selected_real_kind(1) and selected_real_kind(9). Both happen to map
!  onto real*4 / real*8 with gfortran, but the explicit precision/range request is safer across
!  compilers.)
!================================================================================================

module datatype

  implicit none

  integer, parameter, public :: r1 = selected_real_kind(6)
  integer, parameter, public :: r2 = selected_real_kind(15, 307)
  integer, parameter, public :: i8 = selected_int_kind(18)

end module datatype
