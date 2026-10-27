"""Regression test: the reduced example models must reproduce the stored reference values.

Differences between compilers, machines and thread counts are at the level of rounding errors
(about 1e-12); the tolerance of 1e-4 catches every change of the physics or numerics. After an
intended change, regenerate the reference with  python3 tests/regression.py --write.
"""

import json
import os

import numpy as np
import pytest

from conftest import requires_dms
from regression import EXAMPLE_NAMES, REFERENCE, regression_model, summarise

pytestmark = requires_dms

RTOL = 1e-4


@pytest.fixture(scope="module")
def reference():
    if not os.path.isfile(REFERENCE):
        pytest.skip("no reference values (python3 tests/regression.py --write)")
    with open(REFERENCE) as f:
        return json.load(f)


@pytest.mark.parametrize("name", EXAMPLE_NAMES)
def test_example(name, reference, tmp_path):
    ref = reference[name]
    new = summarise(regression_model(name, tmp_path).run())
    assert set(new) == set(ref), "different outputs than in the reference"
    for key, val in ref.items():
        a, b = np.atleast_1d(new[key]), np.atleast_1d(val)
        # quantities that cancel (polarised flux of a face-on disk, U_phi) are compared with an
        # absolute tolerance relative to the scattered light
        scale = max(np.max(np.abs(b)), _scale(ref, key))
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=RTOL * 1e-6 * scale,
                                   err_msg=f"{name}: {key}")


def _scale(ref, key):
    for k in ("img_sca", "sed_scattered") if key in ("img_uphi", "img_qphi", "img_pol",
                                                         "sed_polarised") else ():
        if k in ref:
            return float(np.max(np.abs(ref[k])))
    return 1e-300
