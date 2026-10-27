!________________________________________________________________________________________
!
! 2. 2nd module : physical constants (SI units)
!________________________________________________________________________________________
!
! h, k_B          : CODATA 2018 (exact values in the revised SI)
! c               : exact
! au              : IAU 2012 Resolution B2 (exact)
! pc              : IAU 2015 Resolution B2, 1 pc = (648000/pi) au
! R_sun, M_sun    : IAU 2015 nominal solar radius, M_sun from the nominal GM_sun / G
! G               : CODATA 2018
!________________________________________________________________________________________

module constants

  use datatype

  implicit none

  real(r2), parameter :: pi      = 3.14159265358979323846264338328_r2
  real(r2), parameter :: h       = 6.62607015e-34_r2          ! Planck constant          [J s]
  real(r2), parameter :: k_B     = 1.380649e-23_r2            ! Boltzmann constant       [J/K]
  real(r2), parameter :: c       = 2.99792458e+8_r2           ! speed of light           [m/s]
  real(r2), parameter :: G_grav  = 6.67430e-11_r2             ! gravitational constant   [m^3/kg/s^2]
  real(r2), parameter :: au      = 149597870700.0_r2          ! astronomical unit        [m]
  real(r2), parameter :: pc      = 648000.0_r2 / pi * au      ! parsec                   [m]
  real(r2), parameter :: R_sun   = 6.957e8_r2                 ! nominal solar radius     [m]
  real(r2), parameter :: M_sun   = 1.98840987e30_r2           ! solar mass               [kg]
  real(r2), parameter :: sigma_SB = 5.670374419e-8_r2         ! Stefan-Boltzmann         [W/m^2/K^4]
  real(r2), parameter :: deg2rad = pi / 180.0_r2
  real(r2), parameter :: micron  = 1.0e-6_r2                  ! [m]

end module constants
