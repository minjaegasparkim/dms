# Tests

```bash
make test                         # from the repository root: builds bin/dms and the Python module, runs pytest
python3 -m pytest tests -v        # the same, test by test (after make && make python)
```

| file | content |
|---|---|
| `conftest.py` | shared set-up: a small ring model (`make_model`) and reduced example models |
| `test_physics.py` | checks with known answers: energy balance, image vs. SED, belt masses, size file, polarisation limits, geometry |
| `test_blowout.py` | blow-out threshold (`blowout`, `blowout_release`, `blowout_beta`) |
| `test_tdep.py` | temperature-dependent optical constants (file names, several components, `'nearest'`, `'linear'`, input of version 2.0) |
| `test_dhs.py` | DHS large-grain approximation `dhs_xmax` against the full DHS, spheres, optics cache |
| `test_photosputtering.py` | photosputtering: analytic ice-loss radius and photon flux, reference times, residue, UV file, release orbits, eccentric belts with the ice decided per orbit |
| `test_interface.py` | input and error handling, Python interface, command-line program, example notebooks |
| `test_regression.py` | the eight examples against `reference/examples.json` (relative tolerance 1e-4) |
| `regression.py` | the regression models; `python3 tests/regression.py --write` regenerates the reference |

Tests that compare two runs of the same model use a relative tolerance of 1e-10: the order of
the OpenMP sums changes results by up to ~2e-12 in the far Wien tail on a busy machine.

Regenerate the reference only after a change that alters results on purpose, and explain the
change in `CHANGELOG.md`. Nothing is written into the repository: outputs and optics caches go to
temporary directories.
