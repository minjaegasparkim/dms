# Python interface

The Fortran code is compiled into a Python extension with `numpy.f2py`:

```bash
make python              # needs numpy, gfortran; Python >= 3.12 also needs meson and ninja
pip install -e python    # optional
```

Without `pip install`, add the `python/` directory to `sys.path` or `PYTHONPATH`.

## A first model

[examples/quickstart.py](../examples/quickstart.py) defines a model directly in Python (no input
file) and computes an SED and two images in about ten seconds; the README shows it step by step.
[examples/two_belts_observed.py](../examples/two_belts_observed.py) is a more complete example
(about two minutes): two belts, one of them eccentric, icy grains that sublimate near the star,
blow-out for release at pericentre, a polarised-light image with a PSF and a coronagraph mask,
thermal images with a PSF, grain temperatures and the vertical optical depth.

## Running a model

```python
import dms

m = dms.Model.from_file("examples/debris_ring.nml")      # or dms.Model(star={...}, disk={...}, ...)
m.update(disk={"r_in": 60.0, "eccentricity": 0.1, "omega": 30.0},
         image={"n_pixels": 151, "psf": "gaussian", "psf_fwhm": 0.1})
res = m.run()                    # write_files=False, quiet=True by default
```

Parameter names, units and groups are the same as in the input file (docs/input_reference.md).
Lists are Python lists (`disk={"n_belts": 2, "r_in": [5, 80], ...}`). `m.write("model.nml")`
writes an input file for the command-line program.

## Results

| attribute / method | content |
|---|---|
| `res.wavelength` | observing wavelengths [micron] |
| `res.sed` | dict of arrays [Jy]: `thermal`, `scattered`, `Q`, `U`, `polarised`, `star`, `total` |
| `res.image(name, observed=False, unit="Jy/pixel")` | cube [n_wavelength, ny, nx]; name: `thermal`, `scattered`, `Q`, `U`, `polarised`, `Qphi`, `Uphi`, `poldeg`, `total` |
| `res.images` | names of the available images |
| `res.tau_map`, `res.mass_map_kg` | line-of-sight optical depth and dust mass per pixel |
| `res.offsets("arcsec")` | pixel offsets along one axis |
| `res.f_d`, `res.tau_mean`, `res.tau_peak`, `res.cross_section_au2` | dust proxies |
| `res.tau_radial` | largest radial optical depth in the mid-plane; > 0.1 means the optically thin assumption fails |
| `res.tau_profile` | vertical optical depth: `r_au`, `mean`, `pericentre`, `apocentre` |
| `res.temperature` | `r_au`, `T[n_dist, n_component, n_size]` |
| `res.grains` | `a`, `edges`, `number[n_size, n_component]`, `beta` |
| `res.blowout_beta` | (smallest, largest) blow-out threshold β_crit used; (0, 0) if blow-out is off |
| `res.luminosity_lsun`, `res.pixel_arcsec`, `res.pixel_au` | |
| `res.phase_curve(il, normalise=False)` | (θ [deg], F_sca, F_pol) of the whole disk in 1-degree bins of the scattering angle [Jy per bin] |
| `res.cube("sca")`, `res.cube("reem")`, `res.cube("angle")` | 3-d cubes [il, iz, iy, ix] in Jy, and the scattering angle [iz, iy, ix] (only with `image={"save_3d_cubes": True}`) |

Errors in the input (a missing file, invalid values, ...) raise `dms.DMSError`; the Python
session keeps running. A power-law size grid with fewer than 25 bins per decade gives a
`dms.SizeGridWarning` (input_reference.md, `n_sizes`).

## Speed and parallel runs

- The Mie calculation usually dominates a single run. Set `run={"load_optics": True,
  "save_optics": True, "optics_cache_name": "..."}`: the results are computed once and reused
  for as long as grains, wavelengths and optical data stay the same (e.g. while fitting the
  geometry).
- OpenMP parallelises each run (`OMP_NUM_THREADS`).
- Several processes can share a cache: it is written to a temporary file and renamed, so no
  process reads a partly written file. (`dms.batch` also avoids computing the same cache twice.)
- Screen output of the Fortran code (`quiet=False`) goes to the terminal, not to a Jupyter cell.
- The Fortran state is global, so one model runs at a time per Python process. For many models
  at once (MCMC with emcee, grids), use separate processes (`multiprocessing`, MPI) and set
  `OMP_NUM_THREADS=1` in each.

## Many models: `dms.batch`

`dms.batch` runs input files from a queue folder, several at a time (one process per model):

```python
from dms import batch
base = dms.Model.from_file("examples/debris_ring.nml")
jobs = batch.grid(base, {"disk.inclination": [0, 45, 70], "dust.q": [3.0, 3.5]})
batch.write_queue(jobs, "queue")              # one input file per model
batch.run("queue", "runs", workers=8)         # returns when the queue is empty
rows = batch.collect("runs")                  # list of dicts, one per model
```

or from a terminal: `python -m dms.batch run queue runs -j 8 [--watch]`.

- Input files can be added to `queue/` at any time; with `--watch` the runner keeps picking them
  up (a file `queue/STOP` ends it). Each file is moved to `queue/running/`, then `queue/done/`
  or `queue/failed/`. Several runners may share one queue folder.
- Results: `runs/<name>/` holds the DMS output, `input.nml`, `dms.log`, `meta.json` (labels given
  as a third element of a job: `(name, model, {"star": "A"})`) and `result.json`: wavelengths,
  SED, degree of polarisation against the scattering angle (`dop_theta`, from
  `phase_function.dat`) and, with images, `qphi_over_i` = ΣQ_φ/ΣI_sca.
- `optics_cache_name = 'auto'` (default of `write_queue`) names the Mie cache after the grain radii,
  optical data, wavelengths and angles; models that differ only in geometry, star or `q` share it,
  and the runner does not start a second model that would compute the same cache meanwhile.
- `python -m dms.batch status queue` and `python -m dms.batch requeue queue` (failed → queue).

## Composite grains: `dms.emc`

`dms.emc.mix` runs the effective-medium code emc (Ossenkopf 1991; not part of DMS, give its
executable with `emc=...` or `$DMS_EMC`) without its menus and writes a DMS optical-data file:

```python
from dms import emc
C = emc.Component
emc.mix("bruggeman", [C("astrosilicate_draine03_862", 0.5, "sphere"),
                      C("kim2019/pure_ices/composite_aice.txt", 0.5, "needle")],
        "data/optical_constants/mixtures/as_aice_bg.txt")
```

Theories: `'maxwell-garnett'`, `'bruggeman'`, `'uni-mg'`, `'uni-bg'`; shapes `'host'`,
`'sphere'`, `'needle'`, `'platelet'`, `'coated'` (with `mantle=` index) / `'coating'`, `'shape'`
(with `L=`); `'vacuum'` is emc's built-in vacuum (one component at most). Files are converted
between micron (DMS) and cm (emc); emc runs in a temporary directory, so path lengths do not
matter. The result uses the wavelength grid of the first non-vacuum component (or of component
`grid_from=i`), within the range covered by all components; give the component with the finest
sampling of the features of interest. Note: for a grain made only of core and mantle, `'uni-mg'`
and `'uni-bg'` give the same result (both reduce to the coated-sphere formula).

## Other modules

- `dms.io`: read the output files of the command-line program (`read_image`, `read_sed`,
  `read_profile`, `read_optical_depth`)
- `dms.plot`: quick-look plots of a Result (`image`, `sed`, `temperature`; needs matplotlib)
- `dms.stellar`: build model-atmosphere grids (`build_ck04_grid`, `build_phoenix_grid`,
  `build_blackbody_grid`; also `python -m dms.stellar ...`)
