!____________________________________________________________________________________________________
!
! 8. 8th module :  Mie routines from code miex (Wolf & Voshchinnikov 2004)
!____________________________________________________________________________________________________

! ****************************************************************************************************
! MIEX: MIE SCATTERING CODE FOR LARGE GRAINS
!                                                _____________________________________________________
!                                                Contact information:   swolf@mpia.de (Sebastian Wolf)
! ====================================================================================================
! For a size parameter x = 2 pi a / lambda and a complex refractive index m the routine shexqnn2
! returns
!   * Qext, Qsca, Qabs, Qbk, Qpr, albedo, g
!   * the scattering amplitudes S1, S2 at n_angle2 = 2 n_angle - 1 equidistant angles in [0, pi]
!
! Changes with respect to the previous DMS version
!   * The work array ru was allocated with nterms = 3e8 complex numbers (4.8 GB of address space)
!     at every call. It is now allocated with the number of terms actually required (num).
!   * No SAVEd variables, no module variables: the routine is thread-safe and is called from an
!     OpenMP loop in mie_cals / Qabs_Qsca_S11_S12_cals.
!
! Reference: N.V. Voshchinnikov, "Optics of Cosmic Dust", Astrophys. Space Phys. Rev. 12, 1 (2002)
! ****************************************************************************************************
module mie_routines

  use datatype

  implicit none

  private
  public :: shexqnn2, mie_amplitudes, mie_amplitudes_grid, coated_ab, q_from_ab, dhs_points, &
            dhs_efficiencies, dhs_radii

contains

  ! ==================================================================================================
  ! Ratio of derivative to function for Bessel functions of half order with complex argument,
  ! J'(n)/J(n), by downward recursion from n = num. a = 1/sizepar, ri = complex refractive index.
  ! ==================================================================================================

  pure subroutine aa2(a, ri, num, ru)

    real(kind=r2), intent(in)                   :: a
    complex(kind=r2), intent(in)                :: ri
    integer, intent(in)                         :: num
    complex(kind=r2), dimension(:), intent(out) :: ru

    integer          :: i, i1, j, num1
    complex(kind=r2) :: s, s1

    s       = a / ri
    ru(num) = real(num + 1, kind=r2) * s
    num1    = num - 1

    do j = 1, num1
      i     = num - j
      i1    = i + 1
      s1    = i1 * s
      ru(i) = s1 - 1.0_r2 / (ru(i1) + s1)
    end do

  end subroutine aa2

  !===================================================================================================
  ! shexqnn2: efficiencies and (optionally, doSA = .true.) scattering amplitudes S1, S2 at
  ! n_angle2 = 2*n_angle-1 angles; ier = 1 (size parameter too small) or 2 (too many terms)
  !===================================================================================================

  subroutine shexqnn2(ri, sizepar, Qext, Qsca, Qabs, Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, n_angle, &
                      ab, n_ab)

    complex(kind=r2), intent(in)                :: ri
    real(kind=r2), intent(in)                   :: sizepar
    real(kind=r2), intent(out)                  :: Qext, Qsca, Qabs, Qbk, Qpr, albedo, g
    integer, intent(out)                        :: ier
    complex(kind=r2), dimension(:), intent(out) :: SA1, SA2
    logical, intent(in)                         :: doSA
    integer, intent(in)                         :: n_angle
    complex(kind=r2), allocatable, dimension(:,:), intent(out), optional :: ab   ! (n, 1:2) = a_n, b_n
    integer, intent(out), optional              :: n_ab

    integer, parameter :: nterms = 300000000          ! upper limit for the number of terms
    real(kind=r2), parameter :: eps = 1.0e-20_r2      ! accuracy
    real(kind=r2), parameter :: sizeparmin = 1.0e-6_r2

    integer       :: iterm, num, iu0, iu1, iu2, i_angle2, i_angle
    real(kind=r2) :: r_iterm, factor, pi, ax, besJ0, besJ1, besJ2, besY0, besY1, besY2, b, an, &
                     y, ass, w1, qq, fac, an2, P, T, Si, Co, z
    complex(kind=r2) :: ra0, rb0, ra1, rb1, r, ss, s1, s2, s3, s, rr

    real(kind=r2), dimension(0:1)               :: fact
    real(kind=r2), dimension(:), allocatable    :: mu, fpi, fpi0, fpi1, ftau
    complex(kind=r2), allocatable, dimension(:) :: ru

    ier     = 0
    Qext    = 0.0_r2
    Qsca    = 0.0_r2
    Qabs    = 0.0_r2
    Qbk     = 0.0_r2
    Qpr     = 0.0_r2
    albedo  = 0.0_r2
    g       = 0.0_r2
    fact(0) = 1.0_r2
    fact(1) = 1.0e+250_r2
    factor  = 1.0e+250_r2

    if (doSA) then
      SA1(:) = cmplx(0.0_r2, 0.0_r2, kind=r2)
      SA2(:) = cmplx(0.0_r2, 0.0_r2, kind=r2)
    end if

    if (sizepar <= sizeparmin) then
      ier = 1
      print *, "<!> Error in subroutine shexqnn2: Mie scattering limit exceeded,"
      print *, "    size parameter: ", sizepar
      return
    end if

    pi = 4.0_r2 * atan(1.0_r2)
    ax = 1.0_r2 / sizepar
    b  = 2.0_r2 * ax**2
    ss = (0.0_r2, 0.0_r2)
    s3 = (0.0_r2,-1.0_r2)
    an = 3.0_r2

    ! number of terms for aa2 [Loskutov (1971)]
    y   = sqrt(real(ri * conjg(ri), kind=r2)) * sizepar
    num = int(1.25_r2 * y + 15.5_r2)
    if (y < 1.0_r2) then
      num = int(7.5_r2 * y + 9.0_r2)
    else if ((y > 100.0_r2) .and. (y < 50000.0_r2)) then
      num = int(1.0625_r2 * y + 28.5_r2)
    else if (y >= 50000.0_r2) then
      num = int(1.005_r2 * y + 50.5_r2)
    end if

    if (num > nterms) then
      ier = 2
      print *, "<!> Error in subroutine shexqnn2: number of terms required: ", num
      return
    end if

    allocate(ru(1:num+1))
    if (present(ab)) then
      allocate(ab(num + 1, 2))
      ab = (0.0_r2, 0.0_r2)
      if (present(n_ab)) n_ab = 1
    end if
    if (doSA) then
      allocate(mu(1:n_angle), fpi(1:n_angle), fpi0(1:n_angle), fpi1(1:n_angle), ftau(1:n_angle))
    end if

    call aa2(ax, ri, num, ru)

    ! ------------------------------------------------------------------------------------------
    ! FIRST TERM
    ! ------------------------------------------------------------------------------------------
    iterm = 1

    ass = sqrt(pi / 2.0_r2 * ax)
    w1  = 2.0_r2 / pi * ax
    Si  = sin(sizepar) / sizepar
    Co  = cos(sizepar) / sizepar

    besJ0 =  Si / ass
    besY0 = -Co / ass
    iu0   = 0

    besJ1 = ( Si * ax - Co) / ass
    besY1 = (-Co * ax - Si) / ass
    iu1   = 0
    iu2   = 0

    s   = ru(1) / ri + ax
    s1  = s * besJ1 - besJ0
    s2  = s * besY1 - besY0
    ra0 = s1 / (s1 - s3 * s2)                         ! a_1

    s   = ru(1) * ri + ax
    s1  = s * besJ1 - besJ0
    s2  = s * besY1 - besY0
    rb0 = s1 / (s1 - s3 * s2)                         ! b_1

    if (present(ab)) then
      ab(1, 1) = ra0;  ab(1, 2) = rb0
    end if
    r    = -1.5_r2 * (ra0 - rb0)
    Qext = an * real(ra0 + rb0, kind=r2)
    Qsca = an * real(ra0 * conjg(ra0) + rb0 * conjg(rb0), kind=r2)

    if (doSA) then

      do i_angle = 1, n_angle
        mu(i_angle) = cos((real(i_angle, kind=r2) - 1.0_r2) * (pi/2.0_r2) / real(n_angle-1, kind=r2))
      end do

      fpi0(:) = 0.0_r2
      fpi1(:) = 1.0_r2

      r_iterm = real(iterm, kind=r2)
      fac     = (2.0_r2*r_iterm + 1.0_r2) / (r_iterm * (r_iterm + 1.0_r2))

      do i_angle = 1, n_angle
        i_angle2      = 2 * n_angle - i_angle
        fpi(i_angle)  = fpi1(i_angle)
        ftau(i_angle) = r_iterm * mu(i_angle) * fpi(i_angle) - (r_iterm + 1.0_r2) * fpi0(i_angle)
        P             = (-1.0_r2)**(iterm - 1)
        SA1(i_angle)  = SA1(i_angle) + fac * (ra0*fpi(i_angle)  + rb0*ftau(i_angle))
        T             = (-1.0_r2)**iterm
        SA2(i_angle)  = SA2(i_angle) + fac * (ra0*ftau(i_angle) + rb0*fpi(i_angle))
        if (i_angle /= i_angle2) then
          SA1(i_angle2) = SA1(i_angle2) + fac * (ra0*fpi(i_angle)*P  + rb0*ftau(i_angle)*T)
          SA2(i_angle2) = SA2(i_angle2) + fac * (ra0*ftau(i_angle)*T + rb0*fpi(i_angle)*P)
        end if
      end do

      iterm   = iterm + 1
      r_iterm = real(iterm, kind=r2)

      do i_angle = 1, n_angle
        fpi1(i_angle) = ((2.0_r2*r_iterm - 1.0_r2) / (r_iterm - 1.0_r2)) * mu(i_angle) * fpi(i_angle)
        fpi1(i_angle) = fpi1(i_angle) - r_iterm * fpi0(i_angle) / (r_iterm - 1.0_r2)
        fpi0(i_angle) = fpi(i_angle)
      end do

    else

      iterm = 2

    end if

    ! ------------------------------------------------------------------------------------------
    ! 2., 3., ... num
    ! ------------------------------------------------------------------------------------------
    z = -1.0_r2

    do
      an  = an + 2.0_r2
      an2 = an - 2.0_r2

      if (iu1 == iu0) then
        besY2 = an2 * ax * besY1 - besY0
      else
        besY2 = an2 * ax * besY1 - besY0 / factor
      end if

      if (abs(besY2) > 1.0e+300_r2) then
        besY2 = besY2 / factor
        iu2   = iu1 + 1
      end if

      ! rbrunngraeber 10/14 (suggested by fkirchschlager): avoid besY2*besJ1 overflow
      besJ2 = besY2 / besY1
      besJ2 = w1 / besY1 + besJ2 * besJ1

      r_iterm = real(iterm, kind=r2)

      s   = ru(iterm) / ri + r_iterm * ax
      s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
      s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
      ra1 = s1 / (s1 - s3 * s2)                       ! a_n

      s   = ru(iterm) * ri + r_iterm * ax
      s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
      s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
      rb1 = s1 / (s1 - s3 * s2)                       ! b_n

      if (present(ab)) then
        ab(iterm, 1) = ra1;  ab(iterm, 2) = rb1
        if (present(n_ab)) n_ab = iterm
      end if
      z  = -z
      rr = z * (r_iterm + 0.5_r2) * (ra1 - rb1)
      r  = r + rr
      ss = ss + (r_iterm - 1.0_r2) * (r_iterm + 1.0_r2) / r_iterm * (ra0 * conjg(ra1) &
              + rb0 * conjg(rb1)) + an2 / r_iterm / (r_iterm - 1.0_r2) * (ra0 * conjg(rb0))
      qq   = an * real(ra1 + rb1, kind=r2)
      Qext = Qext + qq
      Qsca = Qsca + an * real(ra1 * conjg(ra1) + rb1 * conjg(rb1), kind=r2)

      if (abs(qq / qext) < eps) exit

      besJ0 = besJ1
      besJ1 = besJ2
      besY0 = besY1
      besY1 = besY2
      iu0   = iu1
      iu1   = iu2
      ra0   = ra1
      rb0   = rb1

      if (doSA) then

        r_iterm = real(iterm, kind=r2)
        fac     = (2.0_r2 * r_iterm + 1.0_r2) / (r_iterm * (r_iterm + 1.0_r2))

        do i_angle = 1, n_angle
          i_angle2      = 2 * n_angle - i_angle
          fpi(i_angle)  = fpi1(i_angle)
          ftau(i_angle) = r_iterm * mu(i_angle) * fpi(i_angle) - (r_iterm + 1.0_r2) * fpi0(i_angle)
          P             = (-1.0_r2)**(iterm - 1)
          SA1(i_angle)  = SA1(i_angle) + fac * (ra0*fpi(i_angle)  + rb0*ftau(i_angle))
          T             = (-1.0_r2)**iterm
          SA2(i_angle)  = SA2(i_angle) + fac * (ra0*ftau(i_angle) + rb0*fpi(i_angle))
          if (i_angle /= i_angle2) then
            SA1(i_angle2) = SA1(i_angle2) + fac * (ra0*fpi(i_angle)*P  + rb0*ftau(i_angle)*T)
            SA2(i_angle2) = SA2(i_angle2) + fac * (ra0*ftau(i_angle)*T + rb0*fpi(i_angle)*P)
          end if
        end do

        iterm   = iterm + 1
        r_iterm = real(iterm, kind=r2)

        do i_angle = 1, n_angle
          fpi1(i_angle) = ((2.0_r2*r_iterm - 1.0_r2) / (r_iterm - 1.0_r2)) * mu(i_angle) * fpi(i_angle)
          fpi1(i_angle) = fpi1(i_angle) - r_iterm * fpi0(i_angle) / (r_iterm - 1.0_r2)
          fpi0(i_angle) = fpi(i_angle)
        end do

      else

        iterm = iterm + 1

      end if

      if (iterm >= num) exit

    end do

    Qext   = b * Qext
    Qsca   = b * Qsca
    Qbk    = 2.0_r2 * b * real(r * conjg(r), kind=r2)
    Qpr    = Qext - 2.0_r2 * b * real(ss, kind=r2)
    Qabs   = Qext - Qsca
    albedo = Qsca / Qext
    g      = (Qext - Qpr) / Qsca

    deallocate(ru)
    if (doSA) deallocate(mu, fpi, fpi0, fpi1, ftau)

  end subroutine shexqnn2

  !===================================================================================================
  ! Scattering amplitudes at arbitrary angles from the Mie coefficients a_n, b_n (n = 1..nab) of
  ! shexqnn2 (same convention: S1 = sum (2n+1)/(n(n+1)) (a_n pi_n + b_n tau_n),
  !                            S2 = sum (2n+1)/(n(n+1)) (a_n tau_n + b_n pi_n)).
  !===================================================================================================

  pure subroutine mie_amplitudes(ab, nab, mu, S1, S2)

    complex(kind=r2), dimension(:,:), intent(in) :: ab
    integer, intent(in)                          :: nab
    real(kind=r2), dimension(:), intent(in)      :: mu
    complex(kind=r2), dimension(:), intent(out)  :: S1, S2

    integer       :: n, k
    real(kind=r2) :: p0, p1, pn, tn, rn, fac

    S1 = (0.0_r2, 0.0_r2);  S2 = (0.0_r2, 0.0_r2)
    do k = 1, size(mu)
      p0 = 0.0_r2;  p1 = 1.0_r2                      ! pi_0, pi_1
      do n = 1, nab
        rn  = real(n, kind=r2)
        pn  = p1
        tn  = rn * mu(k) * pn - (rn + 1.0_r2) * p0
        fac = (2.0_r2 * rn + 1.0_r2) / (rn * (rn + 1.0_r2))
        S1(k) = S1(k) + fac * (ab(n, 1) * pn + ab(n, 2) * tn)
        S2(k) = S2(k) + fac * (ab(n, 1) * tn + ab(n, 2) * pn)
        p1 = ((2.0_r2 * rn + 1.0_r2) * mu(k) * pn - (rn + 1.0_r2) * p0) / rn
        p0 = pn
      end do
    end do

  end subroutine mie_amplitudes

  !===================================================================================================
  ! Scattering amplitudes at the n_angle2 = 2 n_angle - 1 equidistant angles 0 ... pi (the grid of
  ! shexqnn2) from the coefficients a_n, b_n. Each angle theta < 90 deg is evaluated together with
  ! pi - theta (pi_n(-mu) = (-1)^(n-1) pi_n(mu), tau_n(-mu) = (-1)^n tau_n(mu)).
  !===================================================================================================

  pure subroutine mie_amplitudes_grid(ab, nab, n_angle, S1, S2)

    complex(kind=r2), dimension(:,:), intent(in) :: ab
    integer, intent(in)                          :: nab, n_angle
    complex(kind=r2), dimension(:), intent(out)  :: S1, S2          ! (2 n_angle - 1)

    integer          :: n, k, k2
    real(kind=r2)    :: p0, p1, pn, tn, rn, fac, mu, pi, sgn
    complex(kind=r2) :: ap, bt, at, bp

    pi = 4.0_r2 * atan(1.0_r2)
    S1 = (0.0_r2, 0.0_r2);  S2 = (0.0_r2, 0.0_r2)
    do k = 1, n_angle
      k2 = 2 * n_angle - k
      mu = cos(real(k - 1, kind=r2) * (pi / 2.0_r2) / real(n_angle - 1, kind=r2))
      if (k == n_angle) mu = 0.0_r2
      p0 = 0.0_r2;  p1 = 1.0_r2
      sgn = 1.0_r2                                      ! (-1)^(n-1)
      do n = 1, nab
        rn  = real(n, kind=r2)
        pn  = p1
        tn  = rn * mu * pn - (rn + 1.0_r2) * p0
        fac = (2.0_r2 * rn + 1.0_r2) / (rn * (rn + 1.0_r2))
        ap = ab(n, 1) * pn;  bt = ab(n, 2) * tn
        at = ab(n, 1) * tn;  bp = ab(n, 2) * pn
        S1(k) = S1(k) + fac * (ap + bt)
        S2(k) = S2(k) + fac * (at + bp)
        if (k2 /= k) then
          S1(k2) = S1(k2) + fac * sgn * (ap - bt)
          S2(k2) = S2(k2) + fac * sgn * (bp - at)
        end if
        p1 = ((2.0_r2 * rn + 1.0_r2) * mu * pn - (rn + 1.0_r2) * p0) / rn
        p0 = pn
        sgn = -sgn
      end do
    end do

  end subroutine mie_amplitudes_grid

  !===================================================================================================
  ! Coated sphere (core: refractive index m1, size parameter x1; mantle: m2, x2 > x1) in vacuum.
  ! Coefficients a_n, b_n (n = 1 ... nab) in the convention of shexqnn2 (Bohren & Huffman 1983:
  ! xi_n = psi_n - i chi_n, m = n + i k with k >= 0).
  !
  ! Algorithm of Yang (2003, Appl. Opt. 42, 1710) as used by Pena & Pal (2009, Comput. Phys.
  ! Commun. 180, 2348): logarithmic derivatives D1_n = psi_n'/psi_n (downward recurrence),
  ! D3_n = xi_n'/xi_n from D3_n = D1_n + i / (psi_n xi_n) with the product psi_n xi_n by upward
  ! recurrence, and the ratio
  !     Q_n = [psi_n(m2 x1) / xi_n(m2 x1)] / [psi_n(m2 x2) / xi_n(m2 x2)]
  ! by upward recurrence from Q_0 = exp(2i(z2-z1)) (exp(2i z1) - 1) / (exp(2i z2) - 1), which stays
  ! bounded for absorbing mantles (no overflow of psi_n(m2 x)). For m1 = m2 the result is that of
  ! a homogeneous sphere.
  !===================================================================================================

  subroutine coated_ab(m1, m2, x1, x2, ab, nab, ier)

    complex(kind=r2), intent(in)                              :: m1, m2
    real(kind=r2), intent(in)                                 :: x1, x2
    complex(kind=r2), allocatable, dimension(:,:), intent(out) :: ab
    integer, intent(out)                                      :: nab, ier

    complex(kind=r2), parameter :: ci = (0.0_r2, 1.0_r2)
    complex(kind=r2), allocatable, dimension(:) :: d1c, d1a, d1b, d3a, d3b, d1x
    complex(kind=r2) :: zc, za, zb, pxa, pxb, q, ha, hb, g1, g2, ha1, ha2, xin, xin1, rn_x
    real(kind=r2)    :: psi0, psi1, chi0, chi1, chi2, rn
    integer          :: n, nmx

    ier = 0
    nab = 0
    if (.not. (x2 > 0.0_r2) .or. x1 < 0.0_r2 .or. x1 > x2) then
      ier = 1
      return
    end if
    nab = max(2, int(x2 + 4.05_r2 * x2**(1.0_r2 / 3.0_r2) + 2.0_r2))
    ! start of the downward recurrences: well above max(nab, |m x|) (as in shexqnn2; a margin of
    ! only ~16 terms is not enough for large, nearly non-absorbing spheres)
    nmx = int(1.07_r2 * real(max(nab, int(abs(m2 * x2)), int(abs(m2 * x1)), int(abs(m1 * x1))), kind=r2)) + 30

    zc = m1 * x1;  za = m2 * x1;  zb = m2 * x2
    allocate(d1c(0:nmx), d1a(0:nmx), d1b(0:nmx), d1x(0:nmx), d3a(0:nab), d3b(0:nab))
    call d1_down(zc, nmx, d1c)
    call d1_down(za, nmx, d1a)
    call d1_down(zb, nmx, d1b)
    call d1_down(cmplx(x2, 0.0_r2, kind=r2), nmx, d1x)

    ! D3 by upward recurrence of psi_n xi_n
    pxa = 0.5_r2 * (1.0_r2 - exp(2.0_r2 * ci * za));  d3a(0) = ci
    pxb = 0.5_r2 * (1.0_r2 - exp(2.0_r2 * ci * zb));  d3b(0) = ci
    do n = 1, nab
      rn  = real(n, kind=r2)
      pxa = pxa * (rn / za - d1a(n-1)) * (rn / za - d3a(n-1))
      pxb = pxb * (rn / zb - d1b(n-1)) * (rn / zb - d3b(n-1))
      d3a(n) = d1a(n) + ci / pxa
      d3b(n) = d1b(n) + ci / pxb
    end do

    q = exp(2.0_r2 * ci * (zb - za)) * (exp(2.0_r2 * ci * za) - 1.0_r2) / (exp(2.0_r2 * ci * zb) - 1.0_r2)

    ! Riccati-Bessel functions of the outer size parameter (real argument)
    psi0 = cos(x2);  psi1 = sin(x2)                    ! psi_{-1}, psi_0
    chi0 = -sin(x2); chi1 = cos(x2)                    ! chi_{-1}, chi_0

    allocate(ab(nab, 2))
    do n = 1, nab
      rn = real(n, kind=r2)
      ! Q_n
      q = q * ((d3a(n) + rn / za) / (d1a(n) + rn / za)) * ((d1b(n) + rn / zb) / (d3b(n) + rn / zb))
      ! layer 2 (mantle) on the core
      ha1 = d1c(n)
      g1 = m2 * ha1 - m1 * d1a(n)
      g2 = m2 * ha1 - m1 * d3a(n)
      ha = (g2 * d1b(n) - q * g1 * d3b(n)) / (g2 - q * g1)
      g1 = m1 * ha1 - m2 * d1a(n)
      g2 = m1 * ha1 - m2 * d3a(n)
      hb = (g2 * d1b(n) - q * g1 * d3b(n)) / (g2 - q * g1)
      ! psi_n(x2) from psi_{n-1} and D1 (stable for n > x2), chi_n upward
      ha2  = cmplx(psi1, 0.0_r2, kind=r2)                ! psi_{n-1}
      psi0 = psi1
      psi1 = psi1 / real(d1x(n) + rn / x2, kind=r2)       ! psi_n
      chi2 = (2.0_r2 * rn - 1.0_r2) / x2 * chi1 - chi0
      xin1 = cmplx(psi0, -chi1, kind=r2)                  ! xi_{n-1}
      xin  = cmplx(psi1, -chi2, kind=r2)                  ! xi_n
      chi0 = chi1;  chi1 = chi2
      rn_x = cmplx(rn / x2, 0.0_r2, kind=r2)
      ab(n, 1) = ((ha / m2 + rn_x) * psi1 - ha2) / ((ha / m2 + rn_x) * xin - xin1)
      ab(n, 2) = ((m2 * hb + rn_x) * psi1 - ha2) / ((m2 * hb + rn_x) * xin - xin1)
      if (.not. (abs(ab(n, 1)) < huge(1.0_r2) .and. abs(ab(n, 2)) < huge(1.0_r2))) then
        ier = 2
        return
      end if
    end do

    deallocate(d1c, d1a, d1b, d1x, d3a, d3b)

  contains

    pure subroutine d1_down(z, nm, d)
      complex(kind=r2), intent(in)  :: z
      integer, intent(in)           :: nm
      complex(kind=r2), intent(out) :: d(0:nm)
      integer :: k
      d(nm) = (0.0_r2, 0.0_r2)
      do k = nm, 1, -1
        d(k-1) = real(k, kind=r2) / z - 1.0_r2 / (d(k) + real(k, kind=r2) / z)
      end do
    end subroutine d1_down

  end subroutine coated_ab

  !===================================================================================================
  ! Efficiencies from a_n, b_n (Bohren & Huffman 1983, Sects. 4.4 and 4.5); x = outer size parameter
  !===================================================================================================

  pure subroutine q_from_ab(ab, nab, x, Qext, Qsca, Qbk, Qpr, g)

    complex(kind=r2), dimension(:,:), intent(in) :: ab
    integer, intent(in)                          :: nab
    real(kind=r2), intent(in)                    :: x
    real(kind=r2), intent(out)                   :: Qext, Qsca, Qbk, Qpr, g

    integer          :: n
    real(kind=r2)    :: rn, sg, gq
    complex(kind=r2) :: bk

    Qext = 0.0_r2;  Qsca = 0.0_r2;  gq = 0.0_r2;  bk = (0.0_r2, 0.0_r2);  sg = -1.0_r2
    do n = 1, nab
      rn = real(n, kind=r2)
      Qext = Qext + (2.0_r2 * rn + 1.0_r2) * real(ab(n, 1) + ab(n, 2), kind=r2)
      Qsca = Qsca + (2.0_r2 * rn + 1.0_r2) * (abs(ab(n, 1))**2 + abs(ab(n, 2))**2)
      gq = gq + (2.0_r2 * rn + 1.0_r2) / (rn * (rn + 1.0_r2)) * real(ab(n, 1) * conjg(ab(n, 2)), kind=r2)
      if (n < nab) gq = gq + rn * (rn + 2.0_r2) / (rn + 1.0_r2) * &
                        real(ab(n, 1) * conjg(ab(n+1, 1)) + ab(n, 2) * conjg(ab(n+1, 2)), kind=r2)
      bk = bk + (2.0_r2 * rn + 1.0_r2) * sg * (ab(n, 1) - ab(n, 2))
      sg = -sg
    end do
    Qext = 2.0_r2 / x**2 * Qext
    Qsca = 2.0_r2 / x**2 * Qsca
    gq   = 4.0_r2 / x**2 * gq                           ! g * Qsca
    Qbk  = abs(bk)**2 / x**2
    Qpr  = Qext - gq
    g    = 0.0_r2
    if (Qsca > 0.0_r2) g = gq / Qsca

  end subroutine q_from_ab

  !===================================================================================================
  ! Distribution of hollow spheres (Min et al. 2005, A&A 432, 909): vacuum fractions f uniformly
  ! distributed in [0, fmax]. Gauss-Legendre nodes f(j) and weights w(j) (sum w = 1).
  !===================================================================================================

  pure subroutine dhs_points(fmax, nf, f, w)

    real(kind=r2), intent(in)                :: fmax
    integer, intent(in)                      :: nf
    real(kind=r2), dimension(nf), intent(out) :: f, w

    integer       :: i, j, it
    real(kind=r2) :: z, z1, p1, p2, p3, pp, pi

    pi = 4.0_r2 * atan(1.0_r2)
    do i = 1, nf
      z = cos(pi * (real(i, kind=r2) - 0.25_r2) / (real(nf, kind=r2) + 0.5_r2))
      do it = 1, 100
        p1 = 1.0_r2;  p2 = 0.0_r2
        do j = 1, nf
          p3 = p2;  p2 = p1
          p1 = ((2.0_r2 * j - 1.0_r2) * z * p2 - (j - 1.0_r2) * p3) / real(j, kind=r2)
        end do
        pp = nf * (z * p1 - p2) / (z * z - 1.0_r2)
        z1 = z
        z  = z1 - p1 / pp
        if (abs(z - z1) < 1.0e-15_r2) exit
      end do
      f(i) = 0.5_r2 * fmax * (1.0_r2 - z)              ! [-1, 1] -> [0, fmax]
      w(i) = 1.0_r2 / ((1.0_r2 - z * z) * pp * pp)      ! sum = 1 (half of the [-1,1] weights)
    end do

  end subroutine dhs_points

  !===================================================================================================
  ! Hollow sphere of the DHS with vacuum fraction f and the material volume of a solid sphere of size
  ! parameter x_eq: outer size parameter x_out = x_eq (1-f)^(-1/3), core x_in = x_out f^(1/3)
  !===================================================================================================

  pure subroutine dhs_radii(x_eq, f, x_in, x_out)
    real(kind=r2), intent(in)  :: x_eq, f
    real(kind=r2), intent(out) :: x_in, x_out
    x_out = x_eq / (1.0_r2 - f)**(1.0_r2 / 3.0_r2)
    x_in  = x_out * f**(1.0_r2 / 3.0_r2)
  end subroutine dhs_radii

  !===================================================================================================
  ! DHS efficiencies, relative to the geometrical cross section pi a_eq^2 of the volume-equivalent
  ! solid sphere: Q = sum_j w_j C_j / (pi a_eq^2) (cross sections averaged over f). g is the
  ! scattering-weighted mean asymmetry parameter, Qpr = Qext - g Qsca.
  !===================================================================================================

  subroutine dhs_efficiencies(ri, x_eq, fmax, nf, Qext, Qsca, Qabs, Qpr, g, ier)

    complex(kind=r2), intent(in) :: ri
    real(kind=r2), intent(in)    :: x_eq, fmax
    integer, intent(in)          :: nf
    real(kind=r2), intent(out)   :: Qext, Qsca, Qabs, Qpr, g
    integer, intent(out)         :: ier

    real(kind=r2), dimension(nf) :: f, w
    complex(kind=r2), allocatable, dimension(:,:) :: ab
    real(kind=r2) :: x_in, x_out, qe, qs, qb, qp, gg, s
    integer       :: j, nab

    Qext = 0.0_r2;  Qsca = 0.0_r2;  Qabs = 0.0_r2;  Qpr = 0.0_r2;  g = 0.0_r2;  ier = 0
    call dhs_points(fmax, nf, f, w)
    do j = 1, nf
      call dhs_radii(x_eq, f(j), x_in, x_out)
      call coated_ab((1.0_r2, 0.0_r2), ri, x_in, x_out, ab, nab, ier)
      if (ier /= 0) return
      call q_from_ab(ab, nab, x_out, qe, qs, qb, qp, gg)
      s = w(j) * (x_out / x_eq)**2
      Qext = Qext + s * qe
      Qsca = Qsca + s * qs
      Qpr  = Qpr  + s * qp
      deallocate(ab)
    end do
    Qabs = Qext - Qsca
    if (Qsca > 0.0_r2) g = (Qext - Qpr) / Qsca

  end subroutine dhs_efficiencies

end module mie_routines
