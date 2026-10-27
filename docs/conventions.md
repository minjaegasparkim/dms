# Conventions

## Coordinate frames

- **Detector frame** (xd, yd, z): xd and yd span the image plane; the observer is at +z.
- **Sky frame** (x, y, z): the detector frame rotated by the position angle PA:
  x = xd cos PA + yd sin PA, y = −xd sin PA + yd cos PA.
- **Disk frame** (X, Y, Z): Z is the disk normal and Y is the line of nodes, i.e. the projected
  major axis. With inclination i: X = x cos i − z sin i, Y = y, Z = x sin i + z cos i.
  +X points to the far side of the disk and −X to the near side, which faces the observer and is
  where forward scattering appears.

For PA = 0 the major axis lies along the image y axis (as in DMS versions before 2026). PA rotates
the major axis counter-clockwise from +y. If the image is displayed with +y = north and
+x = west (east to the left), this is the astronomical position angle, counted from north
through east.

## Eccentric belts

The dust orbits of a belt are confocal (the star sits at one focus), share the eccentricity e and
are aligned (argument of pericentre ω, measured in the disk plane from +Y towards +X). A point at
distance s from the star and azimuth ψ = atan2(X, Y) lies on the orbit with semi-major axis

    a = s (1 + e cos(ψ − ω)) / (1 − e²).

`r_in` and `r_out` bound a. The dust therefore extends from the pericentre r_in(1 − e) to the
apocentre r_out(1 + e). For such a set of orbits the surface density at fixed a does not depend
on azimuth: the time a particle spends per unit azimuth (∝ s²) cancels the area element
s ds dψ = (s²/a) da dψ. DMS therefore uses

    ρ(s, ψ) = n(a) · a / s,

where n is the radial law of the belt. Grain temperatures follow the true distance s, which
produces pericentre glow.

## Vertical structure

- Wedge: ρ = n(a)·a/s inside |Z| ≤ s sin θ (s = 3-d distance from the star, θ = `opening_angle`).
- Gaussian: ρ = n(a)·a/ϖ · exp(−Z²/2h²) with h = (h/r)·ϖ (ϖ = distance in the disk plane),
  cut at |Z| = `gauss_cut`·h.

## Scattering and polarisation

The scattering angle θ is the angle between the direction from the star to the grain and the
direction from the grain to the observer: cos θ = z/r. The Mie matrix elements are normalised
as S11/(k² C_sca), whose integral over the full sphere is 1.

Stokes Q is positive for polarisation along the image x axis. With φ = atan2(yd, xd), the
contribution of a grain to Q and U is P·(cos 2φ, sin 2φ), where P = S12-weighted intensity
(negative for polarisation perpendicular to the scattering plane).

**Phase function on the angle grid.** DMS interpolates S11 and S12 linearly in θ between the
`n_angles` grid angles. For large grains the forward diffraction peak (width ≈ 1/x rad, x = 2πa/λ)
and the first diffraction rings are much narrower than the grid step. Interpolating point
values there would give far too much light at small scattering angles, i.e. for edge-on disks and
lines of sight close to the star. DMS therefore

1. averages S11 and S12 over the angle bin of each grid angle in the forward region
   (θ ≲ 30/x), computed from the Mie coefficients with at least 4 samples per ring period;
2. normalises ∫ S11 dΩ = 1, with the exact weights of the linear interpolation that the
   integration uses, by rescaling S11 and S12 only in the forward region θ < 6/x (at least the
   first three grid angles), where the energy missed by the grid lies. For small grains and
   Henyey–Greenstein functions this region is the whole sphere;
3. refines volume elements (images) and integration cells (SED) close to the forward direction,
   until they are smaller than ¼ of the grid step seen from the star and smaller than half the
   width of the narrowest belt.

Checked against an exact integration of the Mie phase function (single grain size, ring at
100 au, 5° opening angle, 1.6 µm, x = 1 to 3927, inclinations 60° and 88°): SED within 1.5 %,
images within 0.8 % for rings 2 au or more wide. The largest SED deviation (−1.5 %) is for
x ≈ 40 edge-on, where the peak (width ≈ 1.5°) is close to the grid step. For an extremely narrow
ring (0.1 au wide) seen edge-on the images deviate by up to 6 %. Beyond the forward region the
grid-point values are used. For a single large grain size the rings are sampled only at the grid
angles there; a size distribution smooths them out.

**Distribution of hollow spheres** (`grain_shape = 'dhs'`; Min et al. 2005, A&A 432, 909). Each
grain of radius a stands for hollow spheres with a central vacuum fraction f, uniformly
distributed in 0 ≤ f ≤ f_max. The material volume (and so the mass) is that of the solid sphere of
radius a: outer radius a(1 − f)^(−1/3), core radius a(1 − f)^(−1/3) f^(1/3). Cross sections and
the scattering matrix (S11, S12, weighted with C_sca) are averaged over f with `dhs_nf`
Gauss–Legendre points; efficiencies refer to πa². The coated-sphere coefficients follow Yang
(2003, Appl. Opt. 42, 1710; as in Peña & Pal 2009). Min et al. (2005) used f_max = 1 (the most
extreme case; 0.4 fitted the measured polarisation of quartz grains best); Min et al. (2016,
A&A 585, A13) used f_max = 0.8. Checks: for f_max → 0 the results equal those of Mie spheres
(10⁻⁵); coated spheres agree with scattnlay (Peña & Pal) to 10⁻¹⁰ (Q) and 10⁻⁷ (S11, S12) for
size parameters 0.5–10000. f_max = 1 is allowed: with 10 points the largest node is f = 0.987
(shell thickness 0.4 % of the outer radius), and such thin shells agree with scattnlay to 10⁻¹⁰
(Q_ext, Q_sca) for outer size parameters up to 2 × 10⁴. The DHS Mie step is about 10–20 times
slower than for spheres.

Azimuthal Stokes parameters: Q_φ = −(Q cos 2φ + U sin 2φ), U_φ = Q sin 2φ − U cos 2φ.
Q_φ > 0 means polarisation perpendicular to the direction to the star, the usual pattern of
optically thin scattering; U_φ ≈ 0 for single scattering. Observers often measure φ from north
through east. The sign of Q_φ is the same in that convention, while U_φ may change sign.

## Fractional luminosity and optical depth

- f_d = L_dust/L* = Σ N πa² ⟨Q_abs⟩* / (4π r²), with ⟨Q_abs⟩* = ∫Q_abs I* dλ / ∫I* dλ.
  It is accumulated during the integration. The integrated thermal SED reproduces it to about
  0.1 % (energy balance).
- Geometrical optical depth (Q_ext = 1): `image_tau` along the line of sight, τ⊥ vertically,
  and τ_mean = Σ N πa² / disk area, with area = Σ_belts π(r_out² − r_in²)√(1 − e²).

## Limits of the model

- Optically thin: no extinction of the stellar light and single scattering. DMS reports the
  largest radial optical depth in the mid-plane and warns if it exceeds 0.1.
- Blow-out: a grain released at distance r from a parent body with semi-major axis a is unbound
  for β ≥ r/(2a): 1/2 for circular orbits, between (1 − e)/2 (release at pericentre) and
  (1 + e)/2 (release at apocentre) for eccentric ones. By default DMS uses 1/2
  (`blowout_release = 'circular'`); `'pericentre'` or `'apocentre'` use the eccentricity of each
  belt, and `blowout_beta` sets any threshold. The threshold is applied at every distance, not
  at every azimuth: where belts overlap in distance the smaller threshold applies.
- Diffraction peak: in reality the light diffracted by grains in front of the star stays within
  ≈ 1/x rad of the star and is accompanied by the extinction of the starlight, which DMS does not
  compute. For edge-on disks, the scattered SED therefore contains forward-scattered light that an
  observation would not separate from the star.

## Mass normalisation

`dust_mass` is the mass of the grains that survive sublimation and blow-out. With several belts,
`belt_mass_fraction` fixes how this surviving mass is split between them.
