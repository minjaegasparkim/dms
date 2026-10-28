# Input reference

DMS reads an input file in Fortran namelist format:

```bash
bin/dms model.nml
```

The file has seven groups. All of them are optional and can come in any order; a parameter that
is not given keeps its default. Comments start with `!`. Strings go in quotes, logicals are
written `.true.` / `.false.`, and lists are comma-separated (`lambda = 1.6, 24., 70.`).

Units: lengths in **au**, angles in **degrees**, masses in **M_sun**, stellar radius in **R_sun**,
distance in **pc**, grain radii and wavelengths in **micron**, temperatures in **K**.

If DMS is started without an argument it asks the same things interactively (`bin/dms`, or
`bin/dms < answers.txt`).

---

## &run — what to compute and where to write it

| name | default | meaning |
|---|---|---|
| `mode` | `'image'` | `'image'`; `'sed'` (fast direct integration, no images); `'both'` (images, plus an SED obtained by summing them) |
| `emission` | `'both'` | `'thermal'`, `'scattering'`, `'both'`, `'polarisation'` (all scattering modes also give Q, U, Q_φ, U_φ) |
| `output_format` | `'fits'` | `'ascii'`, `'fits'`, `'both'` (images; the SED is always written as ASCII and, with FITS, also as FITS) |
| `output_name` | `''` | appended to all output file names (`image_sca.<name>.fits`) |
| `name_with_parameters` | `.false.` | also append R_in, α, q, a_min, e, i to the file names |
| `data_dir` | `'data/optical_constants/'` | directory of the optical-data files |
| `output_dir` | `'results/'` | output directory (created if missing) |
| `cache_dir` | `'cache/'` | directory of saved optical properties |
| `load_optics` | `.false.` | load saved Mie results (used only if radii, wavelengths, angles and n, k are identical; otherwise they are recomputed) |
| `save_optics` | `.false.` | save the Mie results |
| `optics_cache_name` | `'dms_optics'` | file name in `cache_dir` (`.miex` is appended) |
| `save_temperature` | `.true.` | write the temperature distribution |
| `save_grain_sizes` | `.true.` | write grain numbers, cross sections and β |
| `screen_output` | `.true.` | print progress to the screen |
| `write_files` | `.true.` | write output files (the Python interface sets this to `.false.` by default) |

## &star

| name | default | meaning |
|---|---|---|
| `spectrum` | `'blackbody'` | `'blackbody'`, `'file'`, `'grid'` |
| `teff` | 5800 | effective temperature (black body, or model grid) |
| `logg` | 4.5 | log g [cgs] (model grid) |
| `radius` | 1.0 | stellar radius [R_sun] |
| `mass` | 1.0 | stellar mass [M_sun], used for β (blow-out) |
| `distance` | 10.0 | distance [pc] |
| `spectrum_file` | | `spectrum = 'file'`: two columns, wavelength [micron] and flux density [Jy] observed at `distance`; the path is relative to the working directory |
| `grid_dir` | | `spectrum = 'grid'`: directory with `grid_index.txt` (see `data/stellar_models/README.md`) |

Beyond the last tabulated wavelength of a spectrum, a Rayleigh–Jeans continuation is used.

## &disk

The disk is made of `n_belts` belts. The parameters marked *per belt* are lists with one value
per belt. Inclination, position angle and the vertical structure are shared by all belts.

| name | default | meaning |
|---|---|---|
| `n_belts` | 1 | number of belts (≤ 10) |
| `density_law` | `'powerlaw'` | *per belt*: `'powerlaw'`, `'ertel'`, `'gaussian'`, `'broken'`, `'user'` (see below) |
| `r_in`, `r_out` | (required) | *per belt*: inner and outer semi-major axis [au] |
| `eccentricity` | 0 | *per belt*: eccentricity of the confocal orbits, 0 ≤ e < 1 |
| `omega` | 0 | *per belt*: argument of pericentre [deg], measured in the disk plane from the line of nodes (+Y) towards the far side (+X) |
| `belt_mass_fraction` | 1 | *per belt*: relative dust mass (normalised over all belts; surviving mass after sublimation/blow-out) |
| `alpha` | 1.0 | *per belt*: exponent (meaning depends on the law) |
| `alpha_out` | −3.0 | *per belt*: outer exponent (`'broken'`) |
| `r_c` | | *per belt*: characteristic radius [au] (`'broken'`) |
| `r_start`, `r_peak` | | *per belt*: [au] (`'ertel'`) |
| `mu`, `sigma` | | *per belt*: centre and width [au] (`'gaussian'`) |
| `user_par` | 0 | 10 parameters passed to `user_density_fn` (`'user'`) |
| `user_density_symmetric` | `.false.` | set to `.true.` if your user density is symmetric under Y → −Y; DMS may then compute only half of the image |
| `vertical` | `'wedge'` | `'wedge'`: constant density within ±`opening_angle`; `'gaussian'`: exp(−Z²/2h²) with h = `h_over_r`·ϖ |
| `opening_angle` | 5.0 | half opening angle of the wedge [deg] |
| `h_over_r` | 0.05 | aspect ratio of the Gaussian profile |
| `gauss_cut` | 4.0 | the Gaussian profile is cut at this many scale heights |
| `dust_mass` | 1e-8 | total dust mass [M_sun] of the grains between a_min and a_max, after sublimation and blow-out |
| `inclination` | 0 | 0 = face-on, 90 = edge-on |
| `position_angle` | 0 | position angle of the projected major axis, counter-clockwise from the +y axis of the image |

Radial laws (a = semi-major axis):

| law | density |
|---|---|
| `'powerlaw'` | n ∝ (a/au)^(−α) |
| `'ertel'` | n ∝ s^α exp[α(1 − s)], s = (a − r_start)/(r_peak − r_start), 0 for a < r_start (Ertel et al. 2012) |
| `'gaussian'` | n ∝ (a/au)^α exp[−(a − μ)²/(2σ²)] |
| `'broken'` | n ∝ [(a/r_c)^(−2α) + (a/r_c)^(−2α_out)]^(−1/2), with α > 0 and α_out < 0 (Augereau et al. 1999) |
| `'user'` | `user_density_fn` in `src/user_density.f90` (recompile after editing it) |

The law describes the volume density for the wedge and the mid-plane density for the Gaussian
profile. In an eccentric belt the density at distance s is n(a)·a/s, where a is the semi-major
axis of the orbit through that point (docs/conventions.md).

## &image

| name | default | meaning |
|---|---|---|
| `n_pixels` | 101 | pixels per axis (odd; the star sits in the central pixel). The field half-width is the largest distance of dust from the star |
| `max_delta_t` | 2.0 | split a volume element if the temperature varies across it by more than this [K] |
| `n_split_max` | 3 | maximum number of splittings of one volume element; 3–4 is usually enough |
| `add_star` | `.false.` | add the stellar flux to the central pixel of the total image |
| `psf` | `'none'` | `'none'`, `'gaussian'` (FWHM given), `'diffraction'` (FWHM = 1.029 λ/D) |
| `psf_fwhm` | | FWHM [arcsec]: one value for all wavelengths, or one value per wavelength |
| `telescope_diameter` | | D [m] for `psf = 'diffraction'` |
| `mask_radius` | 0 | radius of an inner mask [arcsec]; masked pixels are NaN in the `*_obs` images |
| `image_unit` | `'Jy/pixel'` | unit of the written images and profiles: `'Jy/pixel'`, `'mJy/arcsec2'`, `'Jy/arcsec2'` |
| `save_3d_cubes` | `.false.` | also keep and write the 3-d cubes: the flux that every volume element (N × N × N) sends towards the observer, and the scattering angle of every element. Memory: 16 × n_wavelength × N³ bytes (limit 8 GB) |

A PSF or a mask produces extra files `image_*_obs`: the convolved and masked images, with the
polarisation quantities computed from the convolved Q and U.

## &wavelengths

Give either a list or a logarithmic range:

| name | meaning |
|---|---|
| `lambda` | list of observing wavelengths [micron] |
| `lambda_min`, `lambda_max`, `n_lambda` | `n_lambda` logarithmically spaced wavelengths |

All wavelengths must lie inside the wavelength range shared by all optical-data files.

## &dust

| name | default | meaning |
|---|---|---|
| `n_components` | (required) | number of dust materials (≤ 50) |
| `optical_data` | | one file per component (wavelength [micron], n, k), relative to `data_dir`; sub-directories are allowed (`'kim2019/porous_astrosilicate/as50vac50'`). The wavelengths of all files are merged (union within their common range), so the order of the components does not matter |
| `bulk_density` | | one value per component [g/cm³] |
| `abundance` | 1 | relative abundance per component (normalised). A component that only forms as a sublimation residue gets 0 |
| `size_distribution` | `'powerlaw'` | `'single'`, `'powerlaw'`, `'file'` |
| `a_min`, `a_max` | 1, 1000 | grain radii [micron] (`'single'` uses `a_min`) |
| `q` | 3.5 | n(a) ∝ a^(−q) |
| `n_sizes` | 100 | number of logarithmic size bins (one grain size per bin). Use ≥ 30 bins per decade (the default gives 33 for 1–1000 µm) for Mie polarisation images, where coarser grids leave radial stripes in Q_φ, and for (sub)millimetre fluxes: with 10 bins per decade the thermal SED of 1–1000 µm grains ripples with the period of the size grid by up to ~30 % between 100 and 2000 µm (100 bins: < 4 %). Below 25 bins per decade DMS prints a note and Python issues a `dms.SizeGridWarning`. In earlier versions the default was 30 |
| `size_file` | | `'file'`: columns radius [micron] and relative number of grains in the bin (increasing radii; bin edges are the geometric means of neighbouring radii) |
| `blowout` | `.true.` | remove grains blown out by radiation pressure |
| `blowout_release` | `'circular'` | threshold of the removal: `'circular'` β ≥ 1/2; `'pericentre'` β ≥ (1 − e)/2 and `'apocentre'` β ≥ (1 + e)/2 with the eccentricity e of each belt (where belts overlap in distance, the smaller threshold applies) |
| `blowout_beta` | (unset) | a value > 0 sets the threshold everywhere; cannot be combined with `'pericentre'`/`'apocentre'` |
| `grain_shape` | `'sphere'` | `'sphere'`: homogeneous spheres (Mie); `'dhs'`: distribution of hollow spheres (Min et al. 2005), a simple model of non-spherical grains (docs/conventions.md) |
| `dhs_fmax` | 0.8 | `'dhs'`: largest vacuum volume fraction f of the hollow spheres (0 < f_max ≤ 1; with f_max = 1 the largest node is f = 0.987 for 10 points) |
| `dhs_nf` | 10 | `'dhs'`: number of f values (Gauss–Legendre points in [0, f_max]) |
| `sublimation` | `.false.` | switch sublimation on |
| `t_sublimation` | −1 | per component [K]; a value ≤ 0 means the component never sublimates |
| `sublimation_residue` | 0 | per component: 0 means the grain vanishes, n means it becomes component n |
| `tdep_component` | 0 | component that uses temperature-dependent optical data (0 = none) |
| `tdep_temperatures`, `tdep_files` | | temperatures [K] and files of the data sets; every grain uses the set closest to its own temperature |
| `emt` | `.false.` | create a mixture with an effective-medium rule before the run |
| `emt_matrix`, `emt_inclusion` | | input files (in `data_dir`) |
| `emt_volume_fraction` | | volume fraction of the inclusion (0 < f < 1) |
| `emt_rule` | `'maxwell-garnett'` | `'maxwell-garnett'`, `'bruggeman'`, `'lichtenecker'`, `'looyenga'` |
| `emt_output` | | output file written to `data_dir`; use it as a component in `optical_data` |

## &scattering

| name | default | meaning |
|---|---|---|
| `n_angles` | 181 | scattering angles between 0° and 180° (odd) |
| `phase_function` | `'mie'` | `'mie'`, `'hg'`, `'hg2'` |
| `g1`, `g2`, `w` | 0, 0, 1 | `'hg'`: g = `g1`; `'hg2'`: w·HG(g1) + (1 − w)·HG(g2) |
| `p_max` | 0 | maximum degree of polarisation for HG: P(θ) = −p_max sin²θ/(1 + cos²θ) |

With HG phase functions the absorption and scattering cross sections are still taken from Mie
theory.
