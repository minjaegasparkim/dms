!____________________________________________________________________________________________________|
!                                                                                                    |
! 6th module : dynamic sub-sampling of one volume element (SubVol)                                   |
!____________________________________________________________________________________________________|
!====================================================================================================|
! A volume element (cube, centre cx, cy, cz in the detector frame, half size h0) is split            |
! recursively into 8 sub-cubes until, for each sub-cube, one of the following holds:                 |
!   a. the temperature difference across the cube (smallest grain, first component, between the     |
!      minimum and maximum distance of the cube to the star) is below maxdelT, AND the density      |
!      at the 8 corners and the centre is either all non-zero and smooth (|(max+min)/2 - centre|      |
!      <= 1% centre) or all zero (with the cube thicker than the local disk thickness excluded),    |
!   b. the cube lies completely outside the disk (rejected without evaluation),                     |
!   c. the maximum number of splittings Nmax has been reached.                                      |
! The radiation of the final sub-cubes is evaluated at their centres (radiation_point).              |
!                                                                                                    |
! Changes with respect to the previous version                                                       |
!   * The sub-cubes are kept on a stack (depth-first) instead of a queue (at most 7 entries per level).|
!     more than 7*Nmax+1 entries. The old queue was allocated with 1+8^Nmax entries, while a         |
!     complete splitting needs 1+8+...+8^Nmax = (8^(Nmax+1)-1)/7 entries: for Nmax >= 2 the          |
!     arrays overflowed (heap corruption / crash after "got this far SV13").                          |
!   * Exact minimum/maximum distance of a cube to the star (not only its corners), and a wedge       |
!     test, to reject empty cubes without evaluating the density.                                    |
!   * Thin disks: a cube with all 9 samples empty is still split if it is thicker than the disk.     |
!   * Temperatures are read from the pre-computed distance grid (no search through the table).       |
!====================================================================================================|

module volume_splitting

  use datatype
  use constants
  use parameters
  use functions
  use radiation_calculation

  implicit none

contains

  subroutine SubVol(cx, cy, cz, h0, wsym, reem, sca, qq, uu, mass_acc, sig_acc, fd_acc, hist, ph_sca, ph_pol)

    real(r2), intent(in)                  :: cx, cy, cz, h0, wsym
    real(r2), dimension(:), intent(inout) :: reem, sca, qq, uu, hist
    real(r2), intent(inout)               :: mass_acc, sig_acc, fd_acc
    real(r2), dimension(:,:), intent(inout) :: ph_sca, ph_pol

    integer,  parameter :: max_level   = 30          ! absolute limit of splittings
    real(r2), parameter :: f_near_star = 0.05_r2     ! near the star: refine while cube > 0.05 r
    real(r2), dimension(4, 8*(max_level+2)) :: stk   ! x, y, z, h of the cubes on the stack
    integer,  dimension(8*(max_level+2))    :: lev   ! splitting level
    integer  :: top, level, k, kb
    real(r2) :: x, y, z, hh, dmin, dmax, dens(9), dcen, dmx, dmn, T_hi, T_lo, q, ang, thc, min_belt_width
    logical  :: split, near, fwd
    integer  :: nzero
    real(r2), dimension(3, 8), parameter :: sgn = reshape( &
         [ 1.0_r2, 1.0_r2, 1.0_r2,   1.0_r2, 1.0_r2,-1.0_r2,   1.0_r2,-1.0_r2, 1.0_r2, &
           1.0_r2,-1.0_r2,-1.0_r2,  -1.0_r2, 1.0_r2, 1.0_r2,  -1.0_r2, 1.0_r2,-1.0_r2, &
          -1.0_r2,-1.0_r2, 1.0_r2,  -1.0_r2,-1.0_r2,-1.0_r2 ], [3, 8])

    min_belt_width = minval(belt_width(1:n_belts))
    if (.not. (min_belt_width > 0.0_r2)) min_belt_width = huge(1.0_r2)

    top = 1
    stk(:, 1) = [cx, cy, cz, h0]
    lev(1)    = 0

    do while (top > 0)

      x = stk(1, top); y = stk(2, top); z = stk(3, top); hh = stk(4, top)
      level = lev(top)
      top = top - 1

      !-----------------------------------------------------------------------------------------
      ! b. reject cubes outside the disk
      !-----------------------------------------------------------------------------------------
      call box_distances(x, y, z, hh, dmin, dmax)
      if (dmin > r_max_3d .or. dmax < r_in_eff) cycle
      if (cube_outside_wedge(x, y, z, hh)) cycle

      split = .false.

      ! Close to the star a volume element can be large compared with its distance to the star
      ! (e.g. dust down to 0.01 au in an image of +-150 au): there the splitting may continue
      ! beyond Nmax (temperature / boundary criteria only) until the element is smaller than
      ! f_near_star times its distance. Elsewhere Nmax applies as before.
      near = (2.0_r2 * hh > f_near_star * max(dmin, r_in_eff))

      ! Forward scattering: close to the line of sight through the star (scattering angle within
      ! two angle-grid steps) the phase function varies strongly; there the elements are split
      ! (also beyond Nmax) until their angular size seen from the star is below 1/4 grid step
      ! and they are smaller than half the width of the narrowest belt.
      fwd = .false.
      if (do_sca .and. z > 0.0_r2) then
        ang = sqrt(3.0_r2) * hh / max(dmin, r_in_eff)                     ! [rad]
        thc = acos(min(1.0_r2, z / sqrt(x*x + y*y + z*z)))                 ! [rad]
        fwd = (thc < 2.0_r2 * d_angle * deg2rad + ang .and. &
               (ang > 0.25_r2 * d_angle * deg2rad .or. 2.0_r2 * hh > 0.5_r2 * min_belt_width))
      end if

      if (level < Nmax .or. ((near .or. fwd) .and. level < max_level)) then

        do k = 1, 8
          dens(k) = density_fn(x + sgn(1,k)*hh, y + sgn(2,k)*hh, z + sgn(3,k)*hh)
        end do
        dens(9) = density_fn(x, y, z)

        nzero = count(dens <= 0.0_r2)

        if (nzero == 9) then
          ! all samples empty: the disk may still pass between them if the cube is thicker
          ! than the local disk (thin wedge) or wider than a narrow belt it overlaps
          split = (2.0_r2 * hh > sin_open * max(dmin, r_in_eff))
          if (.not. split) then
            do kb = 1, n_belts
              if (dmax >= belt_rin_eff(kb) .and. dmin <= belt_rout_eff(kb) / cos(open_angle) .and. &
                  2.0_r2 * hh > 0.5_r2 * belt_width(kb)) then
                split = .true.
                exit
              end if
            end do
          end if
          if (.not. split) cycle
        else if (nzero > 0) then
          split = .true.                                     ! disk boundary inside the cube
        else
          T_hi = T1_dist(dist_index(dmin))
          T_lo = T1_dist(dist_index(dmax))
          if (abs(T_hi - T_lo) > maxdelT .or. fwd) then
            split = .true.
          else if (level < Nmax) then
            dcen = dens(9)
            dmx  = maxval(dens)
            dmn  = minval(dens)
            if (abs(0.5_r2 * (dmx + dmn) - dcen) > 0.01_r2 * dcen) split = .true.
          end if
        end if

      end if

      if (split) then

        q = 0.5_r2 * hh
        do k = 1, 8
          top = top + 1
          stk(:, top) = [x + sgn(1,k)*q, y + sgn(2,k)*q, z + sgn(3,k)*q, q]
          lev(top)    = level + 1
        end do

      else

        call radiation_point(x, y, z, (2.0_r2 * hh)**3, wsym, reem, sca, qq, uu, &
                             mass_acc, sig_acc, fd_acc, hist, ph_sca, ph_pol)

      end if

    end do

  end subroutine SubVol

end module volume_splitting
