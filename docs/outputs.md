# Output files

All output files are written to `output_dir` (default `results/`). `<ext>` stands for the suffix
built from `output_name` (empty if no name is given). Image files carry the wavelength index:
`_#1`, `_#2`, … in the order of `&wavelengths`.

## Images (`mode = 'image'` or `'both'`)

| file | content | unit |
|---|---|---|
| `image_reem` | thermal re-emission | `image_unit` |
| `image_sca` | scattered stellar light | `image_unit` |
| `image_Q`, `image_U` | Stokes Q, U (Q > 0: polarisation along the image x axis) | `image_unit` |
| `image_pol` | polarised intensity √(Q² + U²) | `image_unit` |
| `image_Qphi`, `image_Uphi` | azimuthal Stokes parameters (docs/conventions.md) | `image_unit` |
| `image_poldeg` | degree of polarisation PI/I (I = dust thermal + scattered, star excluded) | 1 |
| `image_total` | thermal + scattered (+ star in the central pixel if `add_star`) | `image_unit` |
| `image_tau` | line-of-sight geometrical optical depth ∫ n πa² ds per pixel | 1 |
| `image_*_obs` | the same after PSF convolution and masking (only with `psf` or `mask_radius`) | `image_unit` |
| `profile_major_axis<ext>_#i.dat` | cut along the projected major axis of the model images (before PSF convolution; bilinear interpolation): offset [au], offset [arcsec], thermal, scattered, Q, U, PI, Q_φ, total | `image_unit` |

Thermal images are written when `emission` includes thermal emission; the scattered-light and
polarisation images when it includes scattering.

**FITS** (`image_<quantity><ext>.fits`): one cube with NAXIS1 = x, NAXIS2 = y,
NAXIS3 = wavelength (`image_tau` has 2 axes). The header gives the pixel scale (`CDELT1/2` in
arcsec, `PIXAU` in au), the star pixel (`CRPIX1/2`), the wavelengths (`WAVE0001`, …),
`BUNIT`, and the main model parameters, including `FDUST` = f_d and, with blow-out, the
threshold used (`BETACMIN`, `BETACMAX`).

**ASCII** (`image_<quantity><ext>_#i.dat`): header lines starting with `#`, then N rows × N
columns. Row j is the image row iy = j (yd increasing) and column i is ix. With NumPy,
`np.loadtxt(file)[iy-1, ix-1]` gives that pixel.

Pixel (ix, iy) is centred at xd = (ix − (N+1)/2)·Δ, yd = (iy − (N+1)/2)·Δ. The image is shown
as on the sky with north up (+yd) and east left (−xd) when PA is counted from north through east
(docs/conventions.md).

Internally the images are kept in Jy/pixel. Summed over all pixels, they give the SED.

## Flux against scattering angle

`phase_function<ext>.dat` (written whenever scattering is computed, in all modes): the scattered
flux of the whole disk in 1-degree bins of the scattering angle θ (the angle between the
directions star → grain and grain → observer). Column 1 is θ at the bin centre; then, for every
wavelength, F_sca and F_pol in Jy per bin. F_pol > 0 means polarisation perpendicular to the
scattering plane. Summed over all bins, F_sca equals the scattered flux of the SED. Only angles
between about 90° − i and 90° + i (widened by the disk thickness) occur; the other bins are 0.
The curve combines the phase function of the grains with the amount of dust seen at each angle.
In image mode it is computed for every sub-volume; in SED mode it is sampled by the 2-degree
azimuthal grid of the direct integration and is therefore somewhat noisier.

## 3-d cubes (`save_3d_cubes = .true.`)

| file | content |
|---|---|
| `cube_reem<ext>_#i.fits`, `cube_sca<ext>_#i.fits` | flux [Jy] that each volume element (x, y, z; z towards the observer) sends to the observer at wavelength i. Summing along z gives the image in Jy/pixel |
| `cube_scattering_angle<ext>.fits` | scattering angle [deg] at the centre of each volume element |

A histogram of `cube_sca` weighted by `cube_scattering_angle` gives the flux against
scattering angle at the resolution of the volume elements; `phase_function.dat` is the exact
version of that histogram.

## SED

`SED<ext>.dat` (with FITS output also `SED<ext>.fits`): columns
1 wavelength [micron], 2 thermal, 3 scattered, 4 Q, 5 U, 6 polarised flux √(Q² + U²),
7 dust (2 + 3), 8 star, 9 total (7 + 8), all in Jy.

`mode = 'sed'` integrates the disk directly (fast). `mode = 'both'` sums the images.

## Tables

| file | content |
|---|---|
| `optical_depth<ext>.dat` | f_d, total cross section, mean τ, largest radial optical depth in the mid-plane (optically thin check); then the vertical (face-on) optical depth τ⊥ = ∫ n πa² dZ against distance in the disk plane [au]: azimuthal mean, towards pericentre, towards apocentre |
| `Temperature_Distribution<ext>.dat` | distance [au], then T [K] for every component (outer loop) and grain size |
| `grain_size_distribution<ext>.dat` | component, a, bin edges [micron], number of grains N in the whole disk, N πa² [au²], mass [kg], β |
| `summary<ext>.txt` | model parameters, f_d, optical depths, SED, list of files, run time |

## Optical-property cache

`cache/<optics_cache_name>.miex` (binary) holds the Mie results. It is used only when the grain
radii, wavelengths, scattering angles and optical constants all match the current run.
It is written to a temporary file (`<name>.miex.tmp<process id>`) that is then renamed, so
models running in parallel never read a partly written cache.
