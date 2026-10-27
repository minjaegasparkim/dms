!____________________________________________________________________________________________________
!
! User-defined density distribution (density law 5)
!____________________________________________________________________________________________________
!
! Edit user_density_fn below and recompile (make; for Python: make python) to model any
! distribution: azimuthal clumps or arcs, spirals, gaps, asymmetric belts, ...
!
! The function is called for every point of a belt whose density law is 5, after the belt
! limits have been applied (semi-major axis a within [r_in, r_out] of the belt, inside the
! vertical extent). DMS then multiplies the returned value by a/s (eccentric belts, see
! functions.f90) and, for the Gaussian vertical structure, by exp(-Z^2/2h^2).
!
! Arguments
!   a        semi-major axis of the orbit through the point [m] (= distance for circular belts)
!   X, Y, Z  disk-frame coordinates [m]: Z = disk normal, Y = line of nodes (projected major
!            axis of the image for PA = 0), +X = far side of the disk
!   k        belt number
!   par      user_par(1:10) from the input file (&disk user_par = ...)
! The result must be >= 0; only relative values matter (the dust mass fixes the scaling).
!
! Useful quantities: azimuth in the disk plane psi = atan2(X, Y) (0 along +Y, 90 deg along +X).
!
! Note: by default the mirror symmetry of the image (yd -> -yd) is NOT used for law 5. If your
! distribution is symmetric under Y -> -Y, set user_density_symmetric = .true. (&disk) to halve
! the computing time of images.
!
! Example below (default): power law in a, n ~ a^(-par(1)), multiplied by a Gaussian azimuthal
! clump of relative amplitude par(2) at azimuth par(3) [deg] with width par(4) [deg]:
!     rho = (a/au)^(-p1) * [1 + p2 exp(-dpsi^2 / (2 p4^2))]
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
    psi   = atan2(X, Y)                                         ! [rad]
    dpsi  = modulo(psi - par(3) * deg2rad + pi, 2.0_r2 * pi) - pi
    width = max(par(4), 1.0e-3_r2) * deg2rad

    rho = (a / au)**(-par(1)) * (1.0_r2 + par(2) * exp(-0.5_r2 * (dpsi / width)**2))
    ! ---------------------------------------------------------------- to here

    if (Z /= Z .or. k < 0) rho = 0.0_r2                         ! (keeps all arguments referenced)
    if (rho < 0.0_r2) rho = 0.0_r2

  end function user_density_fn

end module user_density
