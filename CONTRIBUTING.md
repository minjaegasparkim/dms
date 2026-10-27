# Contributing to DMS

Bug reports, questions and suggestions are welcome as issues. Please include:

- the input file (namelist) and the DMS version (`CHANGELOG.md`)
- the compiler and its version (`gfortran --version`)
- the complete screen output or error message

## Changing the code

- Free-form Fortran (`.f90`) with `implicit none`; Fortran 2003/2008 features are allowed (the
  `Makefile` compiles with `-std=f2008`). Kinds: `r2` for reals, `i8` for large integers.
- The inner integration (`radiation_calculation.f90`, `volume_splitting.f90`, `functions.f90`)
  runs inside OpenMP loops. These routines must only read global variables.
- Errors: call `dms_fail("message")` and `return`. Do not use `stop`; the Python interface
  must survive invalid input.
- New input parameters: add them to `parameters.f90`, `set_defaults` and `finalize_input`
  (`input_setup.f90`), to the namelist in `read_input_file.f90`, to the interactive questions
  if they are relevant there, and to `docs/input_reference.md`.
- New allocatable arrays in `parameters.f90` must also be released in `dms_reset.f90`.
- Tests: `make test` (or `python3 -m pytest tests`) must pass. New features get a test in
  `tests/`. The regression test compares the examples with `tests/reference/examples.json`; if a
  change alters results on purpose, regenerate the reference with `python3 tests/regression.py
  --write` and explain the change in `CHANGELOG.md`.
- Before a pull request: `make debug` must compile without warnings, and the examples should
  run under both `make debug` and `make`. GitHub Actions runs the same checks on every push.
