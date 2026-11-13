!_______________________________________________________________________________________________|
!                                                                                               |
! Orbit-dependent grain populations (dms 2.2)                                                   |
!_______________________________________________________________________________________________|
!===============================================================================================|
! With photosputtering and ps_orbit = 'release', the ice of a grain depends on its orbit. In an  |
! eccentric belt the dust at one distance from the star belongs to different orbits (pericentre  |
! and apocentre sides), so the population of a point is taken from the orbit through it: each    |
! eccentric belt is divided into ps_orbit_bins bins in semi-major axis a (logarithmic), and a    |
! point of belt k with semi-major axis a takes the population of its bin. Temperatures,          |
! sublimation and blow-out stay functions of the distance r of the point.                         |
!                                                                                               |
! The kernels of a pair (orbit bin ob, distance index id) are stored in "cells" after the n_dist |
! distance entries of j_reem, fl_dist, sig_dist, mfrac_dist and mass_hist:                       |
!     cell index kx = n_dist + ob_off(ob) + id - ob_idlo(ob) + 1,  ob_idlo(ob) <= id <= ob_idhi(ob) |
! Circular belts, and all belts without orbit mode, use the distance entries kx = id as before.  |
!                                                                                               |
! The population of a cell: every grain of initial component ic and size ig ends as component     |
!     dest_lost(ig, ic, id)  if it lost its ice in this orbit bin (ps_lost_ob), else              |
!     dest_keep(ig, ic, id)  (sublimation chain and blow-out at that distance; 0 = removed).      |
!===============================================================================================|

module orbit_populations

  use datatype
  use constants
  use parameters
  use functions

  implicit none

  private
  public :: point_populations, orbit_bin, cell_index, cell_number, rho_sigma

contains

  ! orbit bin of semi-major axis a [m] in belt k (0: belt without orbit bins or a outside)
  pure function orbit_bin(k, a) result(ob)
    integer,  intent(in) :: k
    real(r2), intent(in) :: a
    integer              :: ob, j
    real(r2)             :: t
    ob = 0
    if (ob_count(k) <= 0 .or. a <= 0.0_r2) return
    if (ob_count(k) == 1) then
      ob = ob_first(k)
      return
    end if
    t = log(a / ob_alo(ob_first(k))) / log(ob_ahi(ob_first(k) + ob_count(k) - 1) / ob_alo(ob_first(k)))
    j = int(t * real(ob_count(k), kind=r2)) + 1
    j = max(1, min(ob_count(k), j))
    ob = ob_first(k) + j - 1
  end function orbit_bin


  ! cell of orbit bin ob at distance index id (clamped to the distances the bin can reach)
  pure function cell_index(ob, id) result(kx)
    integer, intent(in) :: ob, id
    integer             :: kx
    kx = n_dist + ob_off(ob) + max(ob_idlo(ob), min(ob_idhi(ob), id)) - ob_idlo(ob) + 1
  end function cell_index


  ! populations at a disk-frame point at distance index id: np entries of (density, kernel index
  ! kx, orbit bin ob); without orbit mode one entry (density_disk, id, 0)
  pure subroutine point_populations(X, Y, Zd, id, np, rho, kx, ob)
    real(r2), intent(in)  :: X, Y, Zd
    integer,  intent(in)  :: id
    integer,  intent(out) :: np
    real(r2), dimension(max_belts + 1), intent(out) :: rho
    integer,  dimension(max_belts + 1), intent(out) :: kx, ob
    integer  :: k, b
    real(r2) :: rk, a, rflat

    if (.not. orbit_mode) then
      np = 1
      rho(1) = density_disk(X, Y, Zd)
      kx(1) = id
      ob(1) = 0
      return
    end if

    np = 0
    rflat = 0.0_r2                    ! belts without orbit bins: one entry at the distance
    do k = 1, n_belts
      if (ob_count(k) > 0) then
        call belt_point(X, Y, Zd, k, rk, a)
        rk = belt_amp(k) * rk
        if (rk <= 0.0_r2) cycle
        b = orbit_bin(k, a)
        np = np + 1
        rho(np) = rk
        ob(np)  = b
        kx(np)  = cell_index(b, id)
      else
        rflat = rflat + belt_amp(k) * density_belt(X, Y, Zd, k)
      end if
    end do
    if (rflat > 0.0_r2) then
      np = np + 1
      rho(np) = rflat
      kx(np)  = id
      ob(np)  = 0
    end if
  end subroutine point_populations


  ! number of grains (size ig, final component ic) per kg of initial dust in kernel entry kx
  ! (distance index id, orbit bin ob; ob = 0: distance entry)
  pure function cell_number(ig, ic, kx, id, ob) result(n)
    integer, intent(in) :: ig, ic, kx, id, ob
    real(r2)            :: n
    integer             :: src, c
    if (ob == 0 .or. kx <= n_dist) then
      n = npm(ig, ic, id)
      return
    end if
    n = 0.0_r2
    do src = 1, n_grain_components
      if (npm_init(ig, src) <= 0.0_r2) cycle
      if (ps_lost_ob(ig, src, ob)) then
        c = dest_lost(ig, src, id)
      else
        c = dest_keep(ig, src, id)
      end if
      if (c == ic) n = n + npm_init(ig, src)
    end do
  end function cell_number


  ! sum over the populations at a disk-frame point of density x cross section per kg (optical depth)
  pure function rho_sigma(X, Y, Zd, id) result(s)
    real(r2), intent(in) :: X, Y, Zd
    integer,  intent(in) :: id
    real(r2)             :: s
    integer  :: np, j
    real(r2), dimension(max_belts + 1) :: rho
    integer,  dimension(max_belts + 1) :: kx, ob
    if (.not. orbit_mode) then
      s = density_disk(X, Y, Zd) * sig_dist(id)
      return
    end if
    call point_populations(X, Y, Zd, id, np, rho, kx, ob)
    s = 0.0_r2
    do j = 1, np
      s = s + rho(j) * sig_dist(kx(j))
    end do
  end function rho_sigma

end module orbit_populations
