!____________________________________________________________________________________________________
!
! User-defined density distribution (density law 5): version used for the cover of the manual
!____________________________________________________________________________________________________
!
! To reproduce the cover, copy this file over src/user_density.f90 and rebuild (make python);
! restore the original afterwards (git checkout src/user_density.f90). See README.md.
!
! par(10) selects the distribution:
!   0  the shipped example of src/user_density.f90, unchanged (azimuthal clump):
!        rho = (a/au)^(-p1) * [1 + p2 exp(-dpsi^2 / (2 p4^2))], clump at azimuth p3 [deg],
!        width p4 [deg]
!   1  belt with a Gaussian gap:
!        rho = (a/au)^(-p1) * [1 - p2 exp(-(a/au - p3)^2 / (2 p4^2))], p3 and p4 in au
! Arguments and conventions as in src/user_density.f90.
!____________________________________________________________________________________________________

module user_density

  use datatype
  use constants

  implicit none

contains

  pure function user_density_fn(a, X, Y, Z, k, par) result(rho)

    real(r2), intent(in)               :: a, X, Y, Z
    integer,  intent(in)               :: k
    real(r2), dimension(10), intent(in) :: par
    real(r2)                           :: rho

    real(r2) :: psi, dpsi, width

    ! ---------------------------------------------------------------- edit from here
    if (nint(par(10)) == 1) then
      rho = (a / au)**(-par(1)) * (1.0_r2 - par(2) * exp(-0.5_r2 * ((a / au - par(3)) / par(4))**2))
    else
      psi   = atan2(X, Y)                                         ! [rad]
      dpsi  = modulo(psi - par(3) * deg2rad + pi, 2.0_r2 * pi) - pi
      width = max(par(4), 1.0e-3_r2) * deg2rad
      rho = (a / au)**(-par(1)) * (1.0_r2 + par(2) * exp(-0.5_r2 * (dpsi / width)**2))
    end if
    ! ---------------------------------------------------------------- to here

    if (Z /= Z .or. k < 0) rho = 0.0_r2                         ! (keeps all arguments referenced)
    if (rho < 0.0_r2) rho = 0.0_r2

  end function user_density_fn

end module user_density
