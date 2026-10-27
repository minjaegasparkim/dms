"""Regression models: the example input files, reduced in size.

`summarise` turns a Result into a few numbers that are stored in tests/reference/examples.json
and compared by tests/test_regression.py. Regenerate the reference values only after an
intended change of the results (and say why in CHANGELOG.md):

    python3 tests/regression.py --write
"""

import argparse
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

REFERENCE = os.path.join(HERE, "reference", "examples.json")

EXAMPLE_NAMES = ["debris_ring", "two_belts", "eccentric_ring", "clump_user_density",
                 "polarisation_hg", "edge_on_disk", "sublimation_ice"]

# images that are compared (summed over the image), if present
IMAGES = ["reem", "sca", "qphi", "uphi", "pol", "total"]


def regression_model(name, tmp):
    from conftest import example_model
    m = example_model(name, tmp)
    # small and fast: fewer size bins, no file output, no optics cache
    m.update(run={"write_files": False, "load_optics": False, "save_optics": False})
    m.update(dust={"n_sizes": min(int(m.params["dust"].get("n_sizes", 100)), 30)})
    return m


def summarise(res):
    out = {"f_d": float(res.f_d), "tau_mean": float(res.tau_mean),
           "wavelength": [float(x) for x in res.wavelength]}
    for k in ("thermal", "scattered", "polarised", "star"):
        out["sed_" + k] = [float(x) for x in res.sed[k]]
    if res.has_image:
        for name in IMAGES:
            if name in res.images:
                out["img_" + name] = [float(x) for x in res.image(name).sum(axis=(1, 2))]
    return out


def compute_all():
    import pathlib
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name in EXAMPLE_NAMES:
            res = regression_model(name, pathlib.Path(tmp)).run()
            results[name] = summarise(res)
            print(f"{name}: f_d = {res.f_d:.6e}", flush=True)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--write", action="store_true", help="write tests/reference/examples.json")
    a = p.parse_args()
    res = compute_all()
    if a.write:
        os.makedirs(os.path.dirname(REFERENCE), exist_ok=True)
        with open(REFERENCE, "w") as f:
            json.dump(res, f, indent=1)
        print("written", REFERENCE)
