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
| `dhs_xmax` | 1000 | `'dhs'`: large grains. The cost of the DHS grows with the size parameter x = 2πa/λ (volume-equivalent radius a). A grain with x > `dhs_xmax` is computed at x = `dhs_xmax` with the imaginary part of the refractive index raised so that k x stays the same (limited to k x = 30, where the grain is opaque); its efficiencies (relative to πa²) are used for the larger grain. This holds in geometrical optics, where absorption depends on the size only through k x. Compared with the full DHS: efficiencies within 1 % for absorbing material; in a debris ring with grains up to 1 mm, temperatures within 0.15 % and fluxes at 33–100 µm within 0.05 % (0.4 % at 10 µm), five times faster (CHANGELOG 2.1.1). Applies to the optical-data grid and, without scattered light (`emission = 'thermal'`), to the observing wavelengths; with scattered light the observing wavelengths use the full DHS (the phase function depends on the size). `0` = full DHS for all sizes (version 2.1.0 and earlier). No effect for `grain_shape = 'sphere'` |
| `sublimation` | `.false.` | switch sublimation on |
| `t_sublimation` | −1 | per component [K]; a value ≤ 0 means the component never sublimates |
| `sublimation_residue` | 0 | per component: 0 means the grain vanishes, n means it becomes component n |
| `photosputtering` | `.false.` | switch ice loss by photosputtering on (version 2.2). A grain of a component with `ps_yield` > 0 loses its ice where the time to sputter it, t_ps = 4 ρ_w a / (m_H₂O Y N_abs), is shorter than the reference time `ps_time` (Grigorieva et al. 2007, A&A 475, 755, their eq. 7 with η = 1); it then becomes a grain of `ps_residue`. ρ_w = `ps_ice_fraction` × `bulk_density` is the ice mass per grain volume, N_abs the absorbed UV photons per unit cross section and time in `ps_band`: ∫ Q(λ) F_λ(r) λ/(hc) dλ with the stellar flux F_λ(r) at the distance r of the grain (plus `ps_isrf` × the band mean of Q). Grains are tested at every size and distance, before sublimation and blow-out; the residue can sublimate afterwards. With blow-out on, grains above the blow-out threshold are not tested (they are removed as they are). The ice-loss radius of every size is in the summary and in `res.photosputtering` |
| `ps_yield` | 0 | per component: sputtering yield Y [H₂O molecules per absorbed photon]; 0 = no photosputtering. Grigorieva et al. (2007) take 10⁻³, between the estimates of Westley et al. (1995, laboratory) and Andersson et al. (2006, molecular dynamics) |
| `ps_ice_fraction` | 1 | per component: mass fraction of water ice (1 for pure ice; for a mixture made with emc, the ice mass fraction of the mixture) |
| `ps_residue` | 0 | per component: 0 means the grain vanishes when its ice is gone, n means it becomes component n (e.g. the porous refractory left behind; give it abundance 0). The residue must not have `ps_yield` > 0 |
| `ps_time` | `'shortest'` | reference time t_ref: `'collision'`: collisional lifetime of Grigorieva et al. (2007, eqs 19–20), t_coll = [(a/s₁)⁻² + (a/s₂)^2.7] / (τΩ) for a < s₂ and (a/s₂)^0.3 / (τΩ) above, s₁ = 1.2 s_½, s₂ = 100 s_½, s_½ the largest radius with β ≥ 1/2 at that distance (the smallest grain if no grain reaches 1/2), τ = `ps_tau`, Ω the Kepler frequency; `'pr'`: Poynting–Robertson drag time r² c / (2β G M*) (their eq. 16; it grows as r² like t_ps, so a size loses its ice everywhere or nowhere); `'age'`: `ps_age`; `'shortest'`: the shortest of the P-R time, the collisional lifetime (if `ps_tau` > 0) and the age (if `ps_age` > 0) |
| `ps_tau` | 0 | vertical geometrical optical depth of the disk for `'collision'` (e.g. `res.tau_peak` of a first run) |
| `ps_age` | 0 | age of the system [yr] for `'age'` |
| `ps_absorption` | `'grain'` | Q(λ) in N_abs: `'grain'` the absorption efficiency of the grain (size and optical data of the component; for a mixture this includes the absorption by the refractory part, an upper limit for the photons absorbed by the ice); `'unity'` Q = 1 |
| `ps_band` | 0.0912, 0.2431 | wavelength band [µm] of the sputtering photons (13.6–5.1 eV, Grigorieva et al. 2007) |
| `ps_uv_file` | | UV spectrum (two columns: λ [µm], F_ν [Jy] at the distance of the star, the format of `&star spectrum_file`) used instead of the stellar model where it covers `ps_band`. Model atmospheres and black bodies lack the chromospheric UV of K and M stars, which can exceed the photospheric UV by orders of magnitude; DMS writes a note for stars below 6000 K without such a file |
| `ps_orbit` | `'local'` | `'local'`: the stellar UV at the distance of the grain. `'release'`: the stellar UV averaged over the Kepler orbit the grain follows after its release at r from its parent orbit, and the P-R time of that orbit. The parent orbit follows `blowout_release` (a_p = r for circular orbits, r/(1−e) at pericentre, r/(1+e) at apocentre, e of the belt); the grain orbit has 1/a_g = (1/a_p − 2β/r)/(1−β) and 1−e_g² = a_p(1−e_p²)/((1−β)a_g), its mean UV is r²⟨r⁻²⟩ = r²/(a_g²√(1−e_g²)) times the local value ((1−2β)^{3/2}/(1−β) for circular parents: 0.80, 0.58, 0.15 at β = 0.1, 0.2, 0.4), and its P-R time a_g² c (1−e_g²)^{3/2} / (β G M* (2+3e_g²)) (Wyatt & Whipple 1950). The collisional lifetime and the age are unchanged; grains stay at their release distance. **Eccentric belts** (e > 0) are divided into `ps_orbit_bins` bins in semi-major axis and the ice is decided per orbit: parent orbit (a, e of the belt), released at its pericentre, at its apocentre or from a circular orbit (`blowout_release`), Kepler frequency at a for the collisional lifetime. A point of the belt takes the population of the orbit through it, so grains of one orbit have the same ice state all round the belt; temperatures, sublimation and blow-out stay functions of the distance. Circular belts keep the decision per distance. `res.photosputtering['a_ice_loss']` gives the largest semi-major axis with ice loss of every size |
| `ps_orbit_bins` | 32 | orbit bins (logarithmic in semi-major axis) per eccentric belt for `ps_orbit = 'release'`. The kernels of every (orbit bin, distance) pair are stored: memory ≈ 8 × (number of wavelengths + 4) bytes per cell; DMS stops above 8 GB |
| `ps_isrf` | 0 | interstellar UV photons in the band [photons cm⁻² s⁻¹] added at all distances |
| `tdep` | `.false.` | per component: temperature-dependent optical data (any number of components). For such a component `optical_data` is a base name, e.g. `'kim2019/pure_ices/composite_cice'`, and DMS uses every file `<base>_<T>K` or `<base>_<T>K.<extension>` in `data_dir` (T integer or decimal, leading zeros allowed: `_10K`, `_010K`, `_12.5K`; any extension or none; a lower-case `k` is accepted too). T is read from the file name and the sets are sorted by it. Two files with the same T (e.g. `_010K.txt` and `_10K.dat`) are an error; no file is an error; a single file is used without temperature dependence (with a warning). The files found are listed on the screen and in the summary; `dms.tdep_files(base, data_dir)` gives the same list in Python. Grains colder or hotter than all sets use the coldest or hottest set; DMS warns once per component, with the temperature range of these grains (screen and summary) |
| `tdep_mode` | `'nearest'` | `'nearest'`: every grain uses the set whose temperature is closest to its own (as in version 2.0); `'linear'`: n and k are interpolated linearly in T between the two sets around the grain temperature (see `tdep_dt`). The choice is made for every grain size, component and distance, also in the temperature calculation, in the radiation pressure (β) and for sublimation |
| `tdep_dt` | 5 | `'linear'`: DMS adds interpolated data sets between neighbouring sets, spaced by at most `tdep_dt` [K], and every grain uses the closest of them, so n, k are interpolated to within `tdep_dt`/2 of the grain temperature. Each set is one more Mie calculation: fewer, wider steps are faster |
| `tdep_component` | 0 | input of version 2.0 (still accepted, one component): component that uses temperature-dependent optical data (0 = none); cannot be combined with `tdep` |
| `tdep_temperatures`, `tdep_files` | | version 2.0: temperatures [K] and files of the data sets of `tdep_component` (any order) |
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
