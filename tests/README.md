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
| `test_interface.py` | input and error handling, Python interface, command-line program |
| `test_regression.py` | the seven examples against `reference/examples.json` (relative tolerance 1e-4) |
| `regression.py` | the regression models; `python3 tests/regression.py --write` regenerates the reference |

Regenerate the reference only after a change that alters results on purpose, and explain the
change in `CHANGELOG.md`. Nothing is written into the repository: outputs and optics caches go to
temporary directories.
