!____________________________________________________________________________________________________
!
! 11th module : tools for observation-like images
!____________________________________________________________________________________________________
!
!   azimuthal_stokes : Q_phi, U_phi from Q, U
!       Q_phi = -(Q cos 2phi + U sin 2phi),  U_phi = Q sin 2phi - U cos 2phi,  phi = atan2(yd, xd)
!       With the DMS convention (Q > 0 along the image x axis) Q_phi > 0 for linear polarisation
!       perpendicular to the direction to the star (centro-symmetric pattern of optically thin
!       dust scattering); U_phi ~ 0 for single scattering. (Observers often measure phi from
!       north through east; the sign of Q_phi is the same, U_phi may change sign.)
!   pol_degree       : P = PI / I (only where I > 0)
!   gauss_convolve   : convolution of every wavelength slice with a normalised, truncated (4 sigma)
!                      2-d Gaussian (separable; flux conserving apart from flux leaving the field)
!____________________________________________________________________________________________________

module obs_tools

  use datatype
  use constants

  implicit none

contains

  subroutine azimuthal_stokes(q, u, qphi, uphi, n)

    integer, intent(in) :: n
    real(r2), dimension(:,:,:), intent(in)  :: q, u          ! (il, ix, iy)
    real(r2), dimension(:,:,:), intent(out) :: qphi, uphi
    integer  :: ix, iy, ic0
    real(r2) :: phi, c2, s2

    ic0 = (n + 1) / 2
    do iy = 1, n
      do ix = 1, n
        if (ix == ic0 .and. iy == ic0) then
          qphi(:, ix, iy) = 0.0_r2
          uphi(:, ix, iy) = 0.0_r2
          cycle
        end if
        phi = atan2(real(iy - ic0, kind=r2), real(ix - ic0, kind=r2))
        c2  = cos(2.0_r2 * phi)
        s2  = sin(2.0_r2 * phi)
        qphi(:, ix, iy) = -(q(:, ix, iy) * c2 + u(:, ix, iy) * s2)
        uphi(:, ix, iy) =   q(:, ix, iy) * s2 - u(:, ix, iy) * c2
      end do
    end do

  end subroutine azimuthal_stokes


  subroutine pol_degree(pi_img, i_img, p)

    real(r2), dimension(:,:,:), intent(in)  :: pi_img, i_img
    real(r2), dimension(:,:,:), intent(out) :: p

    where (i_img > 0.0_r2)
      p = pi_img / i_img
    elsewhere
      p = 0.0_r2
    end where

  end subroutine pol_degree


  subroutine gauss_convolve(img, sigma_pix)

    ! img(il, ix, iy) is convolved in place; sigma_pix(il) in pixels (<= 0.25: no convolution)

    real(r2), dimension(:,:,:), intent(inout) :: img
    real(r2), dimension(:), intent(in)        :: sigma_pix
    integer  :: il, nl, n, m, i, j, k
    real(r2), allocatable :: kern(:), tmp(:,:), src(:,:)

    nl = size(img, 1)
    n  = size(img, 2)
    allocate(tmp(n, n), src(n, n))

    do il = 1, nl
      if (sigma_pix(il) <= 0.25_r2) cycle
      m = max(1, ceiling(4.0_r2 * sigma_pix(il)))
      allocate(kern(-m:m))
      do k = -m, m
        kern(k) = exp(-0.5_r2 * (real(k, kind=r2) / sigma_pix(il))**2)
      end do
      kern = kern / sum(kern)

      src = img(il, :, :)
      ! along x
      tmp = 0.0_r2
      do j = 1, n
        do i = 1, n
          do k = max(-m, 1 - i), min(m, n - i)
            tmp(i, j) = tmp(i, j) + kern(k) * src(i + k, j)
          end do
        end do
      end do
      ! along y
      src = 0.0_r2
      do j = 1, n
        do k = max(-m, 1 - j), min(m, n - j)
          src(:, j) = src(:, j) + kern(k) * tmp(:, j + k)
        end do
      end do
      img(il, :, :) = src
      deallocate(kern)
    end do

    deallocate(tmp, src)

  end subroutine gauss_convolve

end module obs_tools
