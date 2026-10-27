"""Shared set-up for the DMS test suite.

The tests use the Python interface (built with ``make python``) and, for the command-line tests,
the program ``bin/dms`` (built with ``make``). Run them from the repository root:

    make test            # or:  python3 -m pytest tests

All models are small (few pixels, few size bins) so that the whole suite runs in a few minutes.
Nothing is written into the repository: output and optics caches go to temporary directories.
"""

import copy
import os
import sys

import pytest

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "python"))

DATA_DIR = os.path.join(REPO, "data", "optical_constants")
EXAMPLES = os.path.join(REPO, "examples")
DMS_EXE = os.path.join(REPO, "bin", "dms")

try:
    import dms  # noqa: E402
    HAVE_DMS = dms.model._dms is not None
except ImportError:                                         # pragma: no cover
    HAVE_DMS = False

def pytest_configure(config):
    # the test models use few size bins on purpose (speed); do not report this every time
    config.addinivalue_line("filterwarnings", "ignore::dms.SizeGridWarning")


requires_dms = pytest.mark.skipif(not HAVE_DMS, reason="Python extension not built (make python)")
requires_exe = pytest.mark.skipif(not os.path.isfile(DMS_EXE), reason="bin/dms not built (make)")


# a small single-ring model used by many tests
BASE = dict(
    run={"mode": "sed", "emission": "both", "write_files": False, "load_optics": False,
         "save_optics": False, "screen_output": False},
    star={"spectrum": "blackbody", "teff": 5800.0, "radius": 1.0, "mass": 1.0, "distance": 10.0},
    disk={"density_law": "powerlaw", "alpha": 1.0, "r_in": 50.0, "r_out": 100.0,
          "vertical": "wedge", "opening_angle": 5.0, "dust_mass": 1.0e-8,
          "inclination": 60.0, "position_angle": 0.0},
    image={"n_pixels": 41, "n_split_max": 2},
    wavelengths={"lambda": [1.6, 24.0, 70.0, 870.0]},
    dust={"n_components": 1, "optical_data": "astrosilicate_draine03_862", "bulk_density": 3.5,
          "size_distribution": "powerlaw", "a_min": 0.5, "a_max": 1000.0, "q": 3.5, "n_sizes": 20,
          "blowout": True},
    scattering={"n_angles": 181, "phase_function": "mie"},
)


def make_model(**groups):
    """BASE model with some groups updated, e.g. make_model(disk={"inclination": 0})."""
    params = copy.deepcopy(BASE)
    for g, entries in groups.items():
        params.setdefault(g, {}).update(entries)
    params["run"]["data_dir"] = DATA_DIR
    return dms.Model(**params)


def example_model(name, tmp_path, **groups):
    """An example input file, reduced in size, with output and caches in tmp_path."""
    m = dms.Model.from_file(os.path.join(EXAMPLES, name + ".nml"))
    m.update(run={"data_dir": DATA_DIR, "output_dir": str(tmp_path / "results"),
                  "cache_dir": str(tmp_path / "cache"), "screen_output": False})
    m.update(image={"n_pixels": 41, "n_split_max": 2})
    for g, entries in groups.items():
        m.update(**{g: entries})
    return m
