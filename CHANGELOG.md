# Changelog

## 2.2.0 (2026-11) – ice loss by photosputtering

- **New `&dust photosputtering` (default off).** UV photons sputter water molecules off icy grains
  (Grigorieva et al. 2007, A&A 475, 755). A grain of a component with `ps_yield` > 0 loses its ice
  where the sputtering time t_ps = 4 ρ_w a / (m_H₂O Y N_abs) is shorter than a reference time; it
  then becomes a grain of `ps_residue` (or vanishes), as for sublimation. The test is made for
  every grain size and distance, so small grains, which lose their ice first (t_ps ∝ a), are
  treated separately from large ones.
  - N_abs: absorbed photons in 5.1–13.6 eV (`ps_band`) from the stellar spectrum of the run, or from
    a UV spectrum `ps_uv_file` (for the chromospheric UV of K and M stars, which stellar models lack),
    plus an optional interstellar field `ps_isrf`. Q(λ) is the absorption efficiency of the grain
    (`ps_absorption = 'grain'`) or 1.
  - Reference time `ps_time`: the collisional lifetime of Grigorieva et al. (eqs 19–20) with the
    vertical optical depth `ps_tau`, the P-R drag time r² c / (2β G M*), the age `ps_age`, or the
    shortest of these (default).
  - `ps_orbit = 'release'` (default `'local'`): the UV flux averaged over the Kepler orbit a grain
    follows after its release from the parent orbit (circular, or pericentre/apocentre of the belt
    as in `blowout_release`), with the P-R time of that orbit (Wyatt & Whipple 1950). For release
    from circular orbits this lowers the UV by (1−2β)^{3/2}/(1−β): 0.80, 0.58 and 0.15 at β = 0.1,
    0.2 and 0.4. Analytic orbits only; the grains stay at their release distance.
  - **Eccentric belts: ice state per orbit.** With `ps_orbit = 'release'`, every eccentric belt is
    divided into `ps_orbit_bins` bins in semi-major axis (default 32) and the ice is decided per
    orbit (parent orbit a, e of the belt; release as in `blowout_release`; Kepler frequency at a).
    A point of the belt takes the population of the orbit through it, so the grains of one orbit
    keep or lose their ice all round the belt; temperatures, sublimation and blow-out stay
    functions of the distance. The kernels of every (orbit bin, distance) pair are stored after the
    distance entries (new module `src/orbit_populations.f90`); SED, images, belt normalisation,
    optical depth and grain numbers use them. Circular belts are unchanged. Python
    `res.photosputtering['a_ice_loss']`. Checked: with no ice lost the orbit bins reproduce the
    run without photosputtering to 1e-10 (two belts, SED and images); for e → 0 they reproduce
    the circular belt; the ice-loss semi-major axis follows the analytic orbit average; energy
    balance and images against the SED quadrature.
  - Inputs per component: `ps_yield` (Grigorieva et al. take 10⁻³), `ps_ice_fraction` (ice mass
    fraction), `ps_residue`.
  - Output: the largest distance at which each grain size loses its ice (summary; Python
    `res.photosputtering['r_ice_loss']`) and the UV photon flux at 1 au. FITS header card `PHOTOSP`.
  - Checked against the analytic ice-loss radius (age mode, Q = 1) to 1e-3 and against an
    independent integral of the black-body photon flux; tests in `tests/test_photosputtering.py`.
  - Grains above the blow-out threshold are not tested: they leave within an orbit and are removed
    by blow-out as they are.
  - Example 4 (`examples/4_ice_loss.ipynb`): sublimation only, photosputtering against collisions
    and P-R drag (local and release orbits) and against P-R drag alone; the ice-loss radius of every
    grain size and the SEDs; an eccentric belt with the ice decided locally and per orbit.
    Input-file example `examples/photosputtering_ice.nml` (also a regression model).
  - Without photosputtering the results are unchanged (the radiation-pressure β is now computed
    before sublimation; the order of the operations on the grains is otherwise the same).
- Tests: comparisons between two runs of the same model now use a relative tolerance of 1e-10
  (was 1e-12). On a busy machine the order of the OpenMP sums changes results by up to ~2e-12 in
  the far Wien tail (fluxes of 1e-49 Jy), which made the DHS tests fail now and then.
- **Code map** `docs/images/code_map.pdf` (`.svg`, `.png`; source `code_map.dot`): the steps of a run
  with the source files and main routines of each step, the entry points and the Python modules;
  also in the manual (Fig. 4).

## 2.1.1 (2026-11) – faster DHS for large grains

- **New `&dust dhs_xmax` (default 1000), which changes the default results of
  `grain_shape = 'dhs'` slightly.** The DHS cost grows with the size parameter, and the large
  grains at short wavelengths took most of the time (for grains up to 1 mm and optical data from
  0.1 micron, 90 % of it was spent at x > 1000). A grain with x > `dhs_xmax` is now computed at
  x = `dhs_xmax` with k raised so that k x stays the same (limited to k x = 30, opaque); in
  geometrical optics the efficiencies depend on the size only through k x. `dhs_xmax = 0` gives the
  full DHS of version 2.1.0.
  - Efficiencies at x = 3e3–7e4 against the full DHS: Qabs, Qpr, Qsca within 1 % for absorbing
    material (astronomical silicate, forsterite in the UV); for nearly transparent forsterite
    (k ~ 1e-7–1e-5, Qabs < 1e-3) Qabs within ~10 %, the size of the change of the full DHS between
    10 and 40 points in f.
  - Debris ring (A star, 27.75–32.25 au, astronomical silicate + forsterite A/B/C 50–295 K,
    0.264–1000 micron, 30 sizes) against the full DHS: grain temperatures within 0.14 %, fluxes at
    33.5–100 micron within 0.05 %, at 10 micron 0.36 %. With `dhs_xmax = 300`: 0.43 %, 0.36 %,
    1.9 %.
  - Time of that model with 100 sizes on 2 threads: 760 s (before: ~71 min, estimated from the
    30-size runs); with 30 sizes 247 s against 1204 s for the full DHS.
  - Used on the optical-data grid and, without scattered light, at the observing wavelengths. With
    scattered light the observing wavelengths keep the full DHS (the phase function depends on the
    size).
  - Opaque grains above `dhs_xmax` have the same efficiencies at a wavelength; they are computed
    once per wavelength and species.
  - The optics cache records `dhs_xmax` (caches of spheres and of the full DHS remain valid).
- **Coated-sphere routine**: the logarithmic derivatives of real arguments (the outer size
  parameter, the vacuum core of the DHS) are computed in real arithmetic; the DHS is ~15 % faster.
  Only rounding changes.
- `grain_shape = 'sphere'` is unchanged bit for bit (checked with one thread).
- Tests: `tests/test_dhs.py` (DHS against the full DHS for `dhs_xmax` = 1000 and 300, spheres
  unaffected, invalid values, optics cache).
- Documentation: `docs/input_reference.md`, manual (Sect. "Distribution of hollow spheres" and the
  input table), `docs/dhs_xmax_comparison.md` (measurements and comparison tables); the summary
  lists the DHS settings.

## 2.1.0 (2026-11) – temperature-dependent optical constants for several components

- **Any number of components with temperature-dependent optical constants** (`&dust tdep`, one
  logical per component). For such a component `optical_data` is a base name, and DMS uses every
  file `<base>_<T>K` or `<base>_<T>K.<extension>` in `data_dir`: T integer or decimal, leading zeros
  allowed (`_10K`, `_010K`, `_12.5K`), any extension or none, lower-case `k` accepted. The sets are
  sorted by T; two files with the same T are an error, no file is an error, a single file is used
  without temperature dependence (with a warning). Version 2.0 allowed one component, with the
  temperatures and files listed in the input.
- The directory is listed with `opendir`/`readdir` through a small C function
  (`src/dir_listing.c`, called with `ISO_C_BINDING` from the new module `src/tdep_search.f90`);
  the Makefile compiles it with `CC` (default `cc`).
- **`tdep_mode = 'linear'`**: n and k interpolated linearly in T between the two sets around the
  grain temperature. DMS adds interpolated sets spaced by at most `tdep_dt` (default 5 K) and every
  grain uses the closest, so n, k are interpolated to within `tdep_dt`/2. `'nearest'` (default)
  is the behaviour of version 2.0.
- Grains (after sublimation and blow-out) colder or hotter than all sets of their component use the
  coldest or hottest set; DMS now reports their temperature range once per component on the screen
  and in the summary (before, silently).
- The screen output and the summary list the files and temperatures found for every component.
- The input of version 2.0 (`tdep_component`, `tdep_temperatures`, `tdep_files`) is still read and
  gives the same results (the sets may now be given in any order); combining it with `tdep` is an
  error.
- Python: new keys `tdep`, `tdep_mode`, `tdep_dt` (also in the optics-cache key of `dms.batch`,
  together with the files found); `dms.tdep_files(base, data_dir)` lists the sets as DMS finds
  them; `Result.species` gives component, temperature, data file or interpolated, and file of every
  optical species.
- **Memory of the scattering matrix**: the scattering kernels are now computed in place, so that
  two arrays of n_lambda x n_angles x n_sizes x n_species exist at a time instead of four (same
  results). DMS prints the memory needed when it exceeds 1 GB; with many optical species (e.g.
  `'linear'` with a small `tdep_dt`) a run with scattered light can exceed the memory of a laptop,
  and `emission = 'thermal'` does not need these arrays.
- **New optical data**: `data/optical_constants/suto2006/`, crystalline forsterite at 50–295 K for
  the axes A, B, C (Suto et al. 2006), renamed to `fo_suto2006_<axis>_<T>K.dat` for the new file
  names.
- Tests: `tests/test_tdep.py` (file names and edge cases, duplicate temperatures, no or one file,
  the input of version 2.0, several components, the linear mode, temperatures outside the sets).
- Documentation: `docs/input_reference.md`, `docs/python.md`, manual (Sect. "Temperature-dependent
  optical constants").

## 2.0.0 (2026-11) – publication release

### Licence and release
- **Licence.** DMS is distributed under the BSD 3-Clause licence (`LICENSE`). The Mie routine of
  miex (subroutines `aa2` and `shexqnn2` in `src/mie_routines.f90`; Wolf & Voshchinnikov 2004) is
  under the CPC Non-Profit Use Licence and is included with the written permission of S. Wolf;
  `LICENSE`, the module header of `src/mie_routines.f90`, the README and the manual
  ("Citing, licence and acknowledgements") say so. `CITATION.cff` gives the licence, the
  repository and Wolf & Voshchinnikov (2004) as a reference.
- `.zenodo.json` describes the software for the Zenodo archive (title, author, licence, keywords,
  and the related papers Kim et al. 2018 and Wolf & Voshchinnikov 2004).
- **Bulk density of Sil50cice50.** `examples/sublimation_ice.nml` and
  `examples/2_two_belts_observed.ipynb` used 1.95 g cm⁻³ for the 50 % silicate / 50 % ice mixture;
  the volume average is 0.5 × 3.5 + 0.5 × 1.0 = 2.25 g cm⁻³, which both now use. The smallest bound
  icy grain in Example 2 is now 5.1 µm. The reference values of `sublimation_ice` in
  `tests/reference/examples.json` were regenerated.
- The LaTeX source of the manual and the scripts that computed the cover images are no longer in
  the repository; `docs/manual/` holds the PDF and the cover images.

### Getting started and example notebooks
- **Prerequisites and setup check.** `requirements.txt` (NumPy ≥ 1.22, meson, ninja, matplotlib,
  astropy, pytest, JupyterLab); `tools/check_setup.py` and `make check` list the compiler, `make`,
  Python, NumPy, meson/ninja, the optional packages, `bin/dms`, the Python extension for the Python
  that runs the check, and whether `import dms` finds this repository or another installation, and
  print the command that fixes each missing item (standard library only, so it runs before anything
  is installed). README and manual start with the prerequisites for macOS, Linux and Windows (WSL)
  and a table of common build problems.
- **Example notebooks** `examples/0_getting_started.ipynb` (setup check, first SED and image,
  troubleshooting), `1_first_model.ipynb`, `2_two_belts_observed.ipynb` and
  `3_parameter_study.ipynb`, stored with their results so that they can be read on GitHub. They
  replace the scripts `quickstart.py`, `two_belts_observed.py` and `parameter_study.py`. The
  notebooks always import DMS from their own repository, so that an older installation elsewhere
  is not used by mistake. The tests run the code cells of the notebooks in a temporary directory.
- **NumPy 1.x.** The energy-balance test used `np.trapezoid`, which exists only from NumPy 2.0 on;
  it now falls back to `np.trapz`, so the test suite runs with NumPy 1.22 and later, as the
  requirements state.

### Manual cover, Python examples, size bins and apocentre glow
- **Cover of the manual**: eight images computed with DMS 2.0 (scattered and polarised light of a
  narrow ring, an edge-on disk, two belts at three wavelengths, pericentre and apocentre glow of an
  eccentric ring, a belt with a gap and a ring with a clump) in `docs/manual/cover/images/`. The
  manual describes the models in "About the cover"; the gap uses a user-defined density.
- Author name written as Minjae G. Kim (manual, README, `CITATION.cff`).
- **Python quick start** (`examples/quickstart.py`): a model defined in Python without an input
  file; SED and two images in about ten seconds, figure `quickstart.png`. The README starts with it,
  the manual shows it in the Python section, and a test runs it (`test_quickstart_example`; the
  CI now installs matplotlib for this test).
- README rewritten as the front page of the repository: images, the Python example, installation,
  command-line use, table of examples.
- **Second Python example** (`examples/two_belts_observed.py`, about two minutes): a warm inner and
  a cold eccentric outer belt, ice–silicate grains that lose their ice above 110 K, blow-out for
  release at pericentre, Q_φ at 1.6 µm with a PSF and a coronagraph mask, thermal images at 24
  and 870 µm with a PSF, grain temperatures and vertical optical depth; figure
  `docs/images/two_belts_observed.png`. Shown in the README (Example 2), `docs/python.md` and the
  manual. `test_two_belts_observed_example` runs it when `DMS_TEST_SLOW=1` (set in the CI).
- **Third Python example** (`examples/parameter_study.py`, about half a minute on two cores): a
  grid of 12 models (3 grain compositions with 0, 20 and 50 % crystalline ice × G and A star ×
  q = 3.5, 3.0) run with `dms.batch` (`grid` with `linked` parameters, `write_queue`, `run` with
  `workers` and `omp_threads`, `collect`, `dop_at`); a second call computes only missing models.
  Degree of polarisation in the 3 µm ice band against wavelength and scattering angle; figure
  `docs/images/parameter_study.png`. Shown in the README (Example 3), `docs/python.md` and the
  manual; tested by `test_parameter_study_example`.
- **Default number of size bins raised from 30 to 100** (`n_sizes`, `input_setup.f90`). Each size
  bin is represented by one grain size, so with 30 bins for 1–1000 µm (10 per decade) the thermal
  SED ripples with the period of the size grid at 100–2000 µm: up to 29 % from a 400-bin
  calculation for the quick-start model (16 % at 870 µm); 100 bins: 4 %. Runs that relied on the
  default now take about three times longer in the Mie part and are more accurate; set
  `n_sizes = 30` for the previous behaviour. The five examples with `n_sizes = 30` now use 100
  (the regression test still runs them with 30 bins, so the reference values are unchanged).
- **Note on coarse size grids for all emission types.** DMS printed a note below 25 bins per decade
  only for Mie scattering; it now prints it for every run (`Qabs_Qsca_S11_S12_cals.f90`), and the
  Python interface issues `dms.SizeGridWarning` (`check_size_grid` in `model.py`). New tests:
  `test_size_grid_warning`, `test_default_number_of_size_bins`. Documented in
  `docs/input_reference.md`, `docs/python.md` and the manual ("Performance and accuracy").
- **Documentation of apocentre glow corrected** (manual Sect. 2, `examples/eccentric_ring.nml`).
  The text said that eccentric belts show apocentre glow at millimetre wavelengths. This holds for
  the flux integrated across the ring and for images whose beam is wider than the ring: at
  870 micron, `eccentric_ring.nml` gives an apocentre-to-pericentre ratio of 1.27 for the flux per
  unit length. In an image that resolves the ring, the peak surface brightness is higher at
  pericentre at all wavelengths (ratio 0.69 at 870 micron), because the surface density at a given
  semi-major axis does not depend on azimuth. The code was right; only the description changed.

### Blow-out threshold, test suite and CI
- **Blow-out threshold for eccentric belts.** Grains released from a parent body at distance r on an
  orbit with semi-major axis a are unbound for β ≥ r/(2a), i.e. β ≥ 1/2 for circular orbits and
  between (1 − e)/2 (pericentre) and (1 + e)/2 (apocentre) for eccentric ones. New `&dust` inputs:
  `blowout_release = 'circular'` (default, β ≥ 0.5 as before), `'pericentre'` or `'apocentre'`
  (threshold from the eccentricity of each belt; where belts overlap in distance the smaller
  threshold applies), and `blowout_beta` (> 0: any threshold). The threshold is set at every
  distance (`beta_crit_dist`, `particle_abundances.f90`), printed in the summary, written to the
  FITS headers (`BETACMIN`, `BETACMAX`) and returned in Python as `Result.blowout_beta`. With the
  default the results are unchanged: all seven examples agree with the previous version to rounding
  (largest relative difference ~1e-16 from the OpenMP summation order).
- **Test suite** (`tests/`, `make test`, pytest; about 1 minute): physical checks (energy balance,
  image SED against direct quadrature, belt mass fractions, size distribution from a file,
  Rayleigh-limit and Henyey–Greenstein polarisation, position-angle rotation, near-side forward
  scattering, eccentric belt extent, mass scaling), the blow-out threshold, input and error
  handling, the command-line program, and a regression test of the seven examples against
  `tests/reference/examples.json` (relative tolerance 1e-4). `python3 tests/regression.py --write`
  regenerates the reference after an intended change of the results.
- **Continuous integration** (`.github/workflows/tests.yml`): on every push and pull request, a
  debug build runs one example with run-time checks and floating-point traps, then the optimised
  program and Python extension are built and the test suite runs (Ubuntu, Python 3.10 and 3.12).
- `make test` target; `tools/migrate_layout.sh` removed (only needed for the old folder layout).
- Documentation: the code is free-form Fortran (`.f90`) that uses some Fortran 2003/2008 features
  (`error stop`, `open(newunit=)`, `execute_command_line`, `block`); it needs a compiler with
  Fortran 2008 support (README, CONTRIBUTING and the manual said "Fortran 2008" without this
  explanation).

### Fix: wavelength grid of the optical data with several components
- The optical data of all components were interpolated onto the wavelength grid of the first file.
  A finely sampled band of a later component was smoothed: in external mixtures of a refractory
  material (listed first) and water ice, the 3 micron ice band was degraded for coarsely sampled
  refractory data (pyroxene and organics from optool, 4-5 points in 2.5-3.6 micron; up to 20
  percentage points in the polarisation at 3.1 micron), and by about 3 percentage points for
  astrosilicate (31 points). The grid is now the union of the wavelengths of all files within their
  common range (`src/reading_temperature_dependent_data.f90`). Results for one component are
  unchanged. External mixtures computed before must be recomputed; optics caches written before are
  rejected automatically (their wavelength grid differs). Checked against an independent Mie code
  (agreement within 0.5 percentage points for astrosilicate, pyroxene and organics with ice).

### Repository and build
- New layout: `src/` (Fortran), `python/` (Python package `dms`), `examples/`, `data/`, `docs/`,
  `tools/`. Run DMS from the repository root (`bin/dms examples/debris_ring.nml`). The program no
  longer has to be started from `src/` with `../ri-data/` next to it.
- One top-level `Makefile`: `make`, `make debug`, `make serial`, `make python`, `make clean`.
  It builds `build/libdms.a` and the program `bin/dms`. The CMake build was removed.
- Speed comparison with the original 2019 code (same model, `input_image` with Nmax = 1, so that the
  original runs): 2.0 s before, 0.74 s now with 1 thread, 0.26 s with 2 threads. With its own
  setting Nmax = 3 the original code stops with an error (sub-volume buffer overflow); the new code
  needs 0.45 s.
- Removed from the repository: `build/` (compiled files), `results/` (607 old output files),
  `miexres/` (342 cache files in the old format, which the new code cannot read),
  `emtres/`, `.DS_Store` files, the input files `input_image`, `input_SED`, `input_image_emt` (old
  question order), and the `visualization/` scripts (they read the old ASCII outputs; replaced by
  `tools/quicklook.py` and `dms.plot`).
- Optical data: `ri-data/` became `data/optical_constants/` (same file names, so old
  component names still work). `Optical data (Kim et al. 2019)/` became
  `data/optical_constants/kim2019/` with sub-directories without spaces:
  `astrosilicate`, `ice_astrosilicate_MG`, `ice_mantle_astrosilicate_core`,
  `ice_platelets_astrosilicate`, `porous_astrosilicate`, `porous_ices`, `pure_ices`,
  `water_ice_tables` (formerly `ice`). `Readme` became `README.txt`. Astrosilicate files were renamed
  `astrosilicate_draine2003.txt`.
- Not carried over: `ri-data/wave` (identical to `ci50_sil50_BR`), `ri-data/sil10cice90 2`
  (identical to `kim2019/ice_astrosilicate_MG/Sil10cice90`), the duplicate
  `Dust_Porous ices/vacuum.txt`, the `n_k.py` scripts (replaced by
  `tools/plot_optical_constants.py`), and `EMT codes/` (third-party C++ and Fortran
  code, including Ossenkopf's emc; DMS has its own EMT, see `src/dust_mixing.f90`).
- New files: `README.md`, `CHANGELOG.md`, `CONTRIBUTING.md`, `CITATION.cff`, `.gitignore`,
  `docs/` (input reference, outputs, conventions, Python), `data/optical_constants/README.md`,
  `data/stellar_models/README.md`, `tools/migrate_layout.sh`.

### Input
- **Input file (Fortran namelist)**: `bin/dms model.nml`. Groups `&run &star &disk &image
  &wavelengths &dust &scattering`, keyword names, all optional with defaults
  (`src/read_input_file.f90`, `docs/input_reference.md`). The interactive questions remain
  (`bin/dms`), renumbered [1]–[40] for the new options.
- Shared set-up (`src/input_setup.f90`): `set_defaults` and `finalize_input` (unit conversion,
  derived quantities, all checks), used by both input routes.
- Directories can be set in the input (`data_dir`, `output_dir`, `cache_dir`). Spectrum and
  size-distribution files are read relative to the working directory.
- `n_pixels` even → next odd number (as before); `n_split_max` limited to 0–30.

### New physics and geometry
- **Several belts** (`n_belts` ≤ 10): each belt has its own density law, radii, eccentricity,
  argument of pericentre and mass fraction (`belt_normalisation.f90`: belt amplitudes are set
  from the surviving mass of each belt by direct integration).
- **Gaussian vertical structure** (`vertical = 'gaussian'`, `h_over_r`, `gauss_cut`), besides
  the wedge.
- **User-defined density** (`density_law = 'user'`): `src/user_density.f90`, 10 parameters
  `user_par`. The example included is a ring with a Gaussian azimuthal clump.
- **Grain sizes from a file** (`size_distribution = 'file'`).
- **Model-atmosphere grids** (`spectrum = 'grid'`): bilinear interpolation in (T_eff, log g)
  in log F_λ. Grids are built with `python -m dms.stellar` (Castelli & Kurucz 2004, PHOENIX
  Husser et al. 2013, black bodies).

### New outputs
- `image_Qphi`, `image_Uphi` (azimuthal Stokes parameters), `image_poldeg` (PI/I).
- **Simulated observations**: Gaussian PSF (FWHM per wavelength) or diffraction limit
  (1.029 λ/D), inner mask (NaN), `image_*_obs` files.
- Image unit: Jy/pixel, mJy/arcsec² or Jy/arcsec² (`image_unit`).
- Q, U, PI and Q_φ in the radial profiles. The polarisation images are now written whenever
  scattering is computed (in the previous rewrite only for `emission = 'polarisation'`).
- `image_total` is always written.
- **Flux against scattering angle**: `phase_function<ext>.dat` with the scattered and polarised flux
  of the whole disk in 1-degree bins of the scattering angle, accumulated for every sub-volume in
  `radiation_point` (all modes; its sum equals the scattered flux of the SED).
- **3-d cubes** (`&image save_3d_cubes = .true.`, interactive question [19a]): flux of every
  volume element towards the observer (`cube_reem`, `cube_sca`) and its scattering angle
  (`cube_scattering_angle`), as in the `image_*_3D` arrays of the sublimation version.
  Python: `Result.phase_curve`, `Result.cube`.
- The vertical optical-depth profile uses a logarithmic radial grid (400 points), which
  resolves narrow inner belts; it is also kept in memory for Python.

### Examples
- New `examples/eccentric_ring.nml` (e = 0.3 ring, 24 and 870 µm: pericentre glow in the
  mid-infrared, apocentre glow in the millimetre) and `examples/edge_on_disk.nml` (smooth broken
  power law, Gaussian vertical profile, i = 88°, diffraction PSF and mask, polarisation).
- `debris_ring.nml`: 100 size bins instead of 30 (own optics cache name). With 30 bins over
  1–1000 µm, the size-integrated Mie S12 still shows single-size resonances, which appear as
  radial spokes in Q_φ images. DMS now prints a note when Mie scattering is used with fewer than
  25 size bins per decade; docs/input_reference.md recommends ≥ 30 per decade for polarisation.
- `eccentric_ring.nml` has its own optics cache name (`astrosil_eccentric`); it shared one with
  `two_belts.nml` at other wavelengths, so the optics were recomputed whenever the two were run
  alternately.

### Batch runs: cache size limit and recovery
- A Mie cache file holds all optical-data wavelengths, sizes and angles (100–200 MB). In a large
  grid with few repeated grain populations the auto caches filled the disk. `dms.batch.run` now
  keeps the cache folder below `cache_max_gb` (default 10 GB; command line `--cache-max-gb`),
  deleting the least recently used files after each model (never one that a running model
  needs); `cache_max_gb = 0` switches the caches off. Results are the same with and without cache.
- A model whose `result.json` cannot be written (e.g. full disk) now counts as failed.
- `requeue --running` (files left in `queue/running/` after an interruption) and
  `requeue --runs RUNS` (done models without a readable `result.json`).

### Water ice at several temperatures
- `data/optical_constants/mastrapa/temperature/`: amorphous (15–120 K) and crystalline (20–140 K)
  water ice, files as provided, with a README of checks (temperature range of the changes, a step
  in n at 1.29 µm in the crystalline files, comparison with the 50 K files).

### Refractory materials, f_max = 1
- `data/optical_constants/optool/`: amorphous pyroxene Mg0.7Fe0.3SiO3 (Dorschner et al. 1995),
  amorphous carbon (Zubko et al. 1996), organics (Henning & Stognienko 1996) from optool
  (MIT licence), and a pyroxene + carbon mixture as the DIANA dust (0.87 : 0.13 by mass,
  Bruggeman, made with emc). The organics data are cut at 10⁴ µm: at 10⁵ µm, 0.01 µm grains fall
  below the size-parameter limit of the Mie routine.
- `dhs_fmax = 1` is accepted (was: f_max < 1). The Gauss–Legendre nodes lie inside (0, f_max); the
  thinnest shells (f = 0.987) agree with scattnlay to 1e-10.

### Hollow spheres (DHS) and Mastrapa ice data
- `grain_shape = 'dhs'` (`&dust`, with `dhs_fmax`, `dhs_nf`): distribution of hollow spheres
  (Min et al. 2005) at all wavelengths (temperatures, β, scattering). New routines in
  `mie_routines.f90`: `coated_ab` (coated sphere, Yang 2003), `q_from_ab`, `mie_amplitudes_grid`
  (amplitudes on the angle grid, using the symmetry θ ↔ π − θ), `dhs_points`, `dhs_radii`,
  `dhs_efficiencies`; DHS branch in `mie_cals.f90` and `Qabs_Qsca_S11_S12_cals.f90`
  (`dhs_observing`, with the forward-peak treatment for every hollow sphere). Checked against
  homogeneous Mie (miex), scattnlay and the Bohren & Huffman coated-sphere algorithm (PyMieScatt).
  While testing, the starting index of the downward recurrence was raised to 1.07 max(n, |m x|) + 30:
  with only 16 extra terms, large nearly non-absorbing spheres were wrong by up to 0.1 %.
- Optics cache format `DMSMIE05` (grain shape and DHS parameters in the header; older caches are
  recomputed).
- `data/optical_constants/mastrapa/`: water ice at 50 K (amorphous, crystalline; Mastrapa et al.)
  joined with the Kim et al. (2019) data outside the measured range.
- `dms.emc.mix(..., grid_from=i)`: wavelength grid of the result taken from component i.

### Batch runs and composite grains
- `dms.batch` (Python, also `python -m dms.batch`): input files in a queue folder, several models
  in parallel (one process each), `--watch` to keep picking up new files, `STOP` file,
  done/failed folders, `requeue`; per model `result.json` with the SED, the degree of
  polarisation against scattering angle and (images) ΣQ_φ/ΣI; job labels (`meta.json`);
  automatic shared optics cache names, without computing the same cache twice at once.
- `dms.emc`: runs the effective-medium code emc (Ossenkopf 1991, external) without its menus
  (Maxwell-Garnett, Bruggeman, universal MG/BG; spheres, needles, platelets, coated spheres,
  ellipsoids with given L, vacuum) and writes DMS optical-data files. Checked against an
  independent implementation of the formulas: agreement within 1e-4 (emc writes 5 digits).
- `save_optics_cache` (`mie_cals.f90`) writes to a temporary file and renames it (C `rename`
  through ISO_C_BINDING), so parallel runs never read a partly written cache and two runs
  saving the same name do not mix their data.

### Phase function at θ = 0 and edge-on images of narrow rings
Clean build, all examples (command line and Python), all sources compiled as one unit with
strict warnings, documentation links, layout migration. Two problems found and fixed:

- **Value of the phase function at θ = 0 (`phase_function_energy_conservation`).** The whole
  difference between ∫S11_lin dΩ and 1 was put into S11(0). Where the forward peak is resolved,
  this difference is only the quadrature error of the grid (~1e-5), but the first interval has
  a very small solid angle, so S11(0) changed by 30–50 % (e.g. a = 1 µm at 870 µm: 0.125 → 0.085
  sr⁻¹ for a phase function that is flat there). Now S11 and S12 are rescaled only in the forward
  region θ < 6/x (at least three grid angles), which is the whole sphere for small grains and
  Henyey–Greenstein functions. The run summary now prints the largest correction.
- **Edge-on images of very narrow rings (`volume_splitting.f90`).** Near the forward direction
  the volume elements stopped splitting at ¼ grid step (≈ 0.5 au at 100 au), wider than a 0.1-au
  ring, so the ring was sampled at too few points: the edge-on image of that test ring was
  20–22 % too bright for grains of 30–100 µm (the 1.3 % quoted above no longer held for this
  test). Elements in the forward region are now also split while they are wider than half the
  narrowest belt.
- One source line (`read_input_file.f90`) was longer than the 132 characters allowed by the
  standard and compiled only with `-ffree-line-length-none`; it is now split.

Exact semi-analytic check after the fixes (single size, ring at 100 au, 5° opening angle,
1.6 µm, a = 0.3–1000 µm, i.e. x = 1–3927, i = 60° and 88°): SED within 1.5 % (worst: x = 39
edge-on); images within 0.8 % for rings 2 au and 40 au wide; for the 0.1-au ring edge-on the
images are within 6 % (within 0.1 % at i = 60°).

### Blow-out with temperature-dependent data, input checks and packaging
Additional tests: all sources compiled as one unit (`-Wall -Wextra -Wimplicit-procedure
-Wconversion-extra`: no argument mismatches between routines, no single-precision constants);
build and run with gfortran 10.5 and 13.3; eccentric belts against Kepler's laws; PSF convolution
against scipy; resolution study; inclination/position-angle symmetries; temperature-dependent
species selection; failure modes; editable pip install used from another directory.

- **Blow-out with temperature-dependent optical constants (bug).** β was computed with the first
  data set of the component only, so grains were removed (or kept) according to the wrong
  optical constants. Example (component with two data sets, the second one used everywhere): the
  smallest size bin was removed although its β was 0.25, and the 70-µm flux was 11 % too low.
  Now β is computed for every data set, and grains are removed at each distance with the data set
  used there; `beta_rp` reports the largest value over the distances.
- **No grains left.** If blow-out or sublimation removes all grains, DMS now says so ("all grain
  sizes are blown out …") instead of "no dust inside the model space".
- **Dust mass ≤ 0** is rejected (gave empty results).
- **Python package.** `pip install -e python` failed because `pyproject.toml` pointed to the README
  outside the package directory (removed). `Result.image` with an unknown unit gives a clear error;
  `phase_curve(normalise=True)` without scattered light no longer divides by zero.
- **Docs.** README states the tested compilers (gfortran 10.5 and 13.3; the "≥ 9" claim was not
  tested); radial profiles are cuts through the model images before PSF convolution; notes on the
  optics cache with several processes and on screen output in Jupyter.

Checks that passed without changes: mass on the pericentre half of eccentric rings equals Kepler's
time fraction (e = 0.3 / 0.5 / 0.8: 0.3120 / 0.1953 / 0.0525 vs 0.3119 / 0.1955 / 0.0520);
PSF convolution identical to `scipy.ndimage.gaussian_filter` (4e-16; flux leaves the field only
for wide PSFs); unit conversions; inclination 120° and −60° give the 60° image mirrored (1e-7);
PA + 360°, ω + 360° identical; image SED converges to the direct SED (51 → 401 pixels:
0.9987 → 1.0001); temperature-dependent data: the nearest data set is used (1e-12).

### Forward scattering, narrow belts, temperature grid and optically thin check
Tests: 31 input cases (edge cases and invalid inputs) under AddressSanitizer, UndefinedBehaviorSanitizer,
`-fcheck=all` and floating-point traps; independent recalculation with miepython and Python; thread
and repeat consistency. Problems found and fixed:

- **Forward scattering of large grains (physics).** S11 was interpolated between grid angles also
  inside the forward diffraction peak, which is much narrower than the grid step (also the case in
  the 2019 code). Grains seen at small scattering angles (edge-on disks, lines of sight near the star)
  received far too much scattered light: in a test (edge-on, 1.6 µm, grains up to 500 µm) the
  scattered flux was 6.5 times too high (0.102 Jy, now 0.0158 Jy). Now: bin averages of S11, S12
  in the forward region from the Mie coefficients (`mie_amplitudes`, `forward_bins`), an
  energy-conserving value at θ = 0 with normalisation of ∫S11 dΩ (`phase_function_energy_conservation`),
  refinement of volume elements (`volume_splitting.f90`) and SED cells (`computing_SED.f90`) near
  the forward direction. Exact semi-analytic check (edge-on ring, x = 393): images within 1.3 %,
  SED within 0.1 %. Optical-property cache format changed (`DMSMIE04`; old caches are recomputed).
  (The value at θ = 0 and the image accuracy were revised later, see "Phase function at θ = 0
  and edge-on images of narrow rings" above.)
- **Narrow belts missed (sampling).** A belt narrower than a volume element could lie between all
  nine sample points and be skipped; for an edge-on ring 0.1 au wide the central pixels stayed
  empty. Volume elements overlapping a belt are now split while they are wider than half the belt
  width (`belt_width`, `input_setup.f90`).
- **Temperature grid.** The distance grid was linear; for disks spanning several decades in
  distance (e.g. 0.01–150 au: grid step 0.0075 au) the nearest-grid-point temperature near the
  inner edge was off by up to ~17 %. The grid is now logarithmic (constant relative resolution).
- **Optically thin check.** New: largest radial optical depth in the mid-plane (`tau_radial_max`),
  written to `optical_depth.dat` and the summary, with a warning above 0.1. The original
  `input_image` model (dust from 0.01 au, 2.7e-7 M_sun) has τ ≈ 5, so it is not optically thin. The
  `two_belts` example was at τ = 0.10; its dust mass was lowered to 1e-7 M_sun.
- **Python.** NumPy floats were written to the namelist as `np.float64(...)` (Fortran read error,
  e.g. `lambda=np.array([...])`); NumPy integers and booleans are handled as well. The reader
  understands Fortran repeat counts (`2*5.0`). `Result.tau_radial` added.
- **Input.** A size-distribution file with a single row was rejected; an empty file name gave a
  confusing message (now: "file name missing"); an unknown namelist group (e.g. `&disc`) was
  ignored silently (now an error).
- **FITS.** Real header values below 1e-99 lost the `E` of the exponent (invalid card); quotes in
  string values are now doubled as the standard requires.
- **EMT output.** Exponent-safe number format.

Checks that passed without changes: Mie efficiencies and amplitudes vs miepython (≤ 1e-7);
grain temperatures (exact); thermal and scattered fluxes, f_d and cross sections vs an independent
calculation (≤ 0.1 %); β (exact); Stokes Q, U of an inclined, rotated ring vs analytic (≤ 1 %);
1 vs 2 threads (2e-13); 15 consecutive Python runs with an error in between (memory stable).

### Near-star refinement (accuracy fix found in the speed comparison)
- `volume_splitting.f90`: close to the star a volume element can be much larger than its distance to
  the star (e.g. the original `input_image`: dust from 0.01 au in an image of ±150 au, pixel 3 au).
  With the Nmax cap the warm inner dust was strongly under-sampled. Now the splitting continues
  beyond Nmax where the temperature or boundary criterion still asks for it and the element is
  larger than 0.05 × its distance to the star (absolute limit 30 levels). Elsewhere Nmax applies
  as before (no change for the ring examples).
- Test with the original `input_image` model (101 pixels, Nmax = 3), thermal flux at 1 / 10 / 100 µm:
  before 1.1e-30 / 0.20 / 7.44 Jy, now 1.15 / 12.5 / 8.03 Jy; direct SED integration 1.14 / 12.7 /
  8.06 Jy (values after the logarithmic temperature grid, see "Forward scattering, narrow belts,
  temperature grid and optically thin check" above; with the old linear grid the 1-µm value was
  1.37 Jy). Run time 1.1 s.

### Python interface
- `python/dms_api.f90` (wrapped with numpy.f2py), `python/build_extension.py`, package
  `python/dms`: `Model` (from a file or dicts, `update`, `run`, `write`), `Result` (SED, images,
  maps, f_d, τ, temperatures, grain numbers), `DMSError`, `dms.io`, `dms.plot`,
  `dms.stellar`, `dms.namelist`. `python/pyproject.toml` for `pip install -e python`.
- Library mode: errors set a flag and return (`dms_fail` in `functions.f90`) instead of stopping
  the program; every `stop` was replaced. `dms_reset.f90` releases all arrays, so several runs
  can follow each other in one session. Screen output and file output can be switched off.

### Other code changes
- New files: `dms_driver.f90` (the run sequence with an error check after each step),
  `dms_reset.f90`, `input_setup.f90`, `read_input_file.f90`, `user_density.f90`,
  `belt_normalisation.f90`, `obs_tools.f90`, `observation_products.f90`.
- `dms.f90` reads the input-file name from the command line (`-h` for help).
- `functions.f90`: density functions rewritten for belts and vertical profiles
  (`density_belt`, `density_disk`); `dms_fail`, `make_dir` (creates directories only when
  they are missing; paths may contain spaces).
- `computing_images.f90`: file output moved to `observation_products.f90`; field of view and
  rejection tests use the largest 3-d distance of dust (`r_max_3d`, larger than the apocentre for
  the Gaussian profile); mirror symmetry checked for every belt.
- `computing_SED.f90`: 48 latitude cells for the Gaussian profile (16 for the wedge).
- `temperature_distribution_cals.f90`: distance grid up to `r_max_3d`.
- `reading_spectrum.f90`: grid reader; the grid run keeps the requested T_eff.
- `reading_temperature_dependent_data.f90`: wavelength-range check with 1e-9 tolerance (a log grid
  ending exactly at the last data point was rejected because of rounding).
- `summary.f90`: belts, vertical structure, simulated observations; summary file only when
  files are written.
- Character variables for file names: 512 characters (was 250).

### Verification of this release (with `astrosilicate_draine03_862`)
- Energy balance, two eccentric belts, Gaussian profile: integrated thermal SED / L* =
  1.0006 × f_d.
- SED summed from images vs direct integration (two belts, Gaussian): within 0.07 %.
- Belt mass fractions 5 % / 95 % recovered as 4.97 % / 95.03 % from the τ⊥ profile.
- Size distribution from a file identical to the same power law (f_d and SED to 1e-6).
- Model grid of black bodies: on a grid node identical to the black-body run; between nodes
  L differs by 0.7 % (black bodies are not linear in log F_λ between 5500 K and 6000 K).
- Flux against scattering angle: sum over the angle bins = scattered SED flux (to 1e-7); 3-d
  cube summed along the line of sight = scattered image (to 1e-16); angles seen at i = 60°:
  26°–154°; image and SED mode agree.
- PSF convolution conserves the flux (0.088750 vs 0.088755 Jy at 70 µm with a 0.3″ mask).
- Interactive questions, input file and Python interface give consistent results; errors in
  Python raise `DMSError` without ending the session.

### Notes on the optical data (found while reorganising)
- `kim2019/pure_ices/composite_cice_085K.txt` is byte-identical to `composite_cice_075K.txt`.
- `kim2019/ice_astrosilicate_MG/Sil30aice70` is byte-identical to `Sil25aice75`.
- `ri-data/sil100` is identical to `astrosilicate_draine03_237`; `water_aice` and
  `water_cice_55K` are identical to `composite_aice.txt` and `composite_cice_055K.txt`.
These were kept as they are.

## 1.9 (2026-10) – rewrite of the 2019 version

### 1. Bugs fixed in the previous version

| File | Problem | Fix |
|---|---|---|
| volume_splitting.f90 | Sub-volume queue allocated with `1 + 8^Nmax` entries, but a full splitting needs `(8^(Nmax+1)−1)/7`; overflow for Nmax ≥ 2 | Depth-first stack, never more than `7·Nmax+1` entries (allocated `8·(Nmax+2)`) |
| radiation_calculation.f90 | Scattering angle `acos(y/r sin i − z/r cos i)` applied the inclination a second time to detector-frame coordinates | `cos θ = z / r` (observer at +z, star at the origin) |
| temperature lookup | Temperature-table index could be 0 / undefined at the grid edges | Pre-computed T(r) on a distance grid; index clamped to 1…n_dist |
| computing_images / SED | `total_mass`, `Numparsum` used without initialisation | All accumulators zeroed per pixel / per run |
| computing_images | Outer rim cut (model box smaller than the disk) | Box half width = apocentre `Rmax(1+e)`, exact cube–sphere distance tests |
| SED-only mode | `stellar_radiation` not allocated | Stellar intensity is a function (`stellar_radiance`) available in all modes |
| reading_spectrum.f90 | Spectrum deallocated right after reading; `Planck_fn` returned 0 for T* = 0, so a measured spectrum gave no heating | Spectrum kept; intensity interpolated log–log, Rayleigh–Jeans tail; L* and an equivalent T_eff computed |
| mie_cals.f90 | Double allocation in the pre-calculated branch; loops to `n_grain_radius+1` | Single allocation, correct bounds; binary cache with consistency check |
| Qabs_Qsca_S11_S12_cals.f90 | `Qabs_075K` / `Qabs_085K` arrays mixed up | Replaced by the generic temperature-dependent species scheme (no hard-coded arrays) |
| mie_routines.f90 | Work array `ru` allocated with 3·10⁸ complex numbers at every call | Allocated with the number of terms actually needed; no SAVEd state (thread-safe) |
| pre_cals.f90 | Riemann sum `a^−q Δa` whose last bin reached beyond a_max | Log bins, exact ∫a^−q da per bin, geometric bin centre |
| file names | Numbers written with leading blanks into file names | `adjustl`/`trim` everywhere |

### 2. New features

1. **Sublimation** (from the newer version): T_subl per component; a sublimating grain vanishes or becomes another component (e.g. ice/silicate → porous silicate). Components are processed in order of increasing T_subl so a residue can sublimate later. Sublimation-zone radius is printed. The dust mass refers to the *surviving* mass.
2. **Radiation-pressure blow-out**: β = 3πR*²∫Q_pr I* dλ / (4GcM*ρa); grains with β ≥ 0.5 removed (optional); β written to the grain-size file.
3. **FITS output** (`fits_io.f90`, no external library): one cube per quantity (x, y, wavelength), BITPIX −64, pixel scale in arcsec and au, wavelengths `WAVE0001…`, model parameters in the header. Choice: ASCII, FITS or both (question [36]). SED also as FITS.
4. **Henyey–Greenstein phase functions**: single HG (g) or two-component HG, p = w·HG(g1) + (1−w)·HG(g2), with Rayleigh-like polarisation P = −p_max sin²θ/(1+cos²θ). Absorption/scattering cross sections stay those of Mie theory.
5. **Position angle** (disk rotation in the image plane, counter-clockwise from +y; PA = 0 reproduces the old orientation).
6. **Geometrical optical depth** τ = ∫ n πa² ds: line-of-sight map (`image_tau`), vertical τ⊥(r) profile (azimuthal mean, towards pericentre and apocentre), total cross section and mean τ.
7. **Fractional luminosity** f_d = L_d/L* = Σ N πa² ⟨Q_abs⟩* / (4πr²), computed during the integration (checked: equals the integrated thermal SED / L* to 0.1 %).
8. **Eccentric disks**: input = eccentricity e (and argument of pericentre ω if e > 0); R_in, R_out are the inner/outer semi-major axes. Confocal orbits with common ω; density ρ(r,f) = ρ₀(a)·a/r with a = r(1+e cos f)/(1−e²); dust from R_in(1−e) to R_out(1+e). Temperatures follow the true distance (pericentre glow).
9. **Smooth broken power law** density (law 4): n ∝ [(a/R_c)^(−2α_in) + (a/R_c)^(−2α_out)]^(−1/2).
10. **Temperature-dependent optical constants** generalised: any component, any number of data sets (temperature + file); a grain uses the set closest to its own temperature (also inside the temperature calculation). Replaces the ~15 000 hard-coded lines.
11. **Radial profiles** along the projected major axis for every wavelength.
12. **Summary file** `../results/summary<ext>.txt` and wall-clock timing.

### 3. Speed

- OpenMP in the image loop (pixels, dynamic schedule), the SED quadrature, all Mie loops, the temperature table and the kernels.
- Mie: only efficiencies on the optical-data grid (no amplitudes); amplitudes only at the observing wavelengths, with n,k interpolated there (more accurate than interpolating Q and S11).
- Temperatures: r(T) tabulated once (matrix product), inverted onto a 20 000-point distance grid; O(1) lookup instead of searches in 0.5 K steps.
- Emission and scattering kernels per kg of dust pre-computed per distance; when the grain population does not depend on distance a single angle-dependent kernel is used.
- Mirror symmetry (half the image) when PA = 0/180° and (e = 0 or ω = ±90°).
- Empty cubes rejected by exact distance and wedge tests; SED mode integrates the disk directly in (r, latitude, azimuth) instead of a full 3-D box.
- Example (2 threads, 101² pixels, Nmax 3, 4 wavelengths, 30 sizes): integration 0.75 s; total 3–10 s, mostly Mie (cached afterwards).

### 4. Behaviour changes to be aware of

- Constants: CODATA 2018 h, k_B, G; IAU au, pc = 648000/π au, nominal R_sun; M_sun 1.98840987×10³⁰ kg.
- Images are in **Jy/pixel**; Stokes Q > 0 along the image x axis; `image_pol` = polarised intensity √(Q²+U²) (positive).
- Image box half width = R_out(1+e) (was R_out).
- Image files: one file per wavelength `_#<index>`; ASCII images are N rows (y) × N columns (x).
- Optical-data files are read until end of file (headers/comment lines skipped); no line count needed.
- `switch_parameter_dependent` is replaced by questions [26] (temperature-dependent data) and [27] (sublimation).
- EMT writes the mixture to `../ri-data/<name>` so it can be used as a normal component.
- Optical-property cache: `../miexres/<name>.miex` (binary); used only if radii, wavelengths, angles and n,k are identical, otherwise recomputed with a warning.
- Single grain size: question [29a] gives the radius; no distribution.

