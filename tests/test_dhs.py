"""Distribution of hollow spheres (grain_shape = 'dhs') and the large-grain approximation
dhs_xmax (DMS 2.1.1)."""

import numpy as np
import pytest

from conftest import make_model, requires_dms

pytestmark = requires_dms

import dms  # noqa: E402


def dhs_model(**dust):
    """Grains up to 100 micron: on the optical-data grid (from 0.05 micron) x reaches ~1e4."""
    m = make_model(dust={"grain_shape": "dhs", "a_min": 1.0, "a_max": 100.0, "n_sizes": 8, **dust})
    m.update(wavelengths={"lambda": [2.0, 24.0, 70.0, 160.0]})
    return m


def test_large_grain_approximation_matches_full_dhs():
    """dhs_xmax = 1000 (default) against the full DHS (dhs_xmax = 0): temperatures and fluxes agree
    to < 0.2 %; with dhs_xmax = 300 to < 1 %."""
    full = dhs_model(dhs_xmax=0.0).run()
    for xmax, tol in ((1000.0, 2e-3), (300.0, 1e-2)):
        res = dhs_model(dhs_xmax=xmax).run()
        np.testing.assert_allclose(res.temperature["T"], full.temperature["T"], rtol=tol)
        np.testing.assert_allclose(res.sed["thermal"], full.sed["thermal"], rtol=tol)
        np.testing.assert_allclose(res.sed["scattered"], full.sed["scattered"], rtol=tol)
        assert res.f_d == pytest.approx(full.f_d, rel=tol)
    default = dhs_model().run()
    np.testing.assert_allclose(default.sed["thermal"], dhs_model(dhs_xmax=1000.0).run().sed["thermal"],
                               rtol=1e-10)


def test_sphere_not_affected():
    """dhs_xmax has no effect on homogeneous spheres. (Two runs differ only by the order of the
    OpenMP sums: ~1e-14, up to ~2e-12 in the far Wien tail when the machine is busy.)"""
    a = make_model(dust={"n_sizes": 8}).run()
    b = make_model(dust={"n_sizes": 8, "dhs_xmax": 300.0}).run()
    np.testing.assert_allclose(a.sed["thermal"], b.sed["thermal"], rtol=1e-10)
    np.testing.assert_allclose(a.sed["scattered"], b.sed["scattered"], rtol=1e-10)
    np.testing.assert_allclose(a.temperature["T"], b.temperature["T"], rtol=1e-10)


@pytest.mark.parametrize("xmax", [-1.0, 5.0])
def test_invalid_xmax(xmax):
    with pytest.raises(dms.DMSError, match="dhs_xmax"):
        dhs_model(dhs_xmax=xmax).run()


def test_cache_records_xmax(tmp_path):
    """The optics cache stores dhs_xmax: a cache made with another value is not used."""
    cache = {"cache_dir": str(tmp_path) + "/", "optics_cache_name": "dhs"}
    m = dhs_model(dhs_xmax=1000.0)
    m.update(run={**cache, "save_optics": True})
    first = m.run()
    m = dhs_model(dhs_xmax=1000.0)
    m.update(run={**cache, "load_optics": True})
    np.testing.assert_allclose(m.run().sed["thermal"], first.sed["thermal"], rtol=1e-10)
    other = dhs_model(dhs_xmax=300.0)
    other.update(run={**cache, "load_optics": True})
    direct = dhs_model(dhs_xmax=300.0).run()
    np.testing.assert_allclose(other.run().sed["thermal"], direct.sed["thermal"], rtol=1e-10)
    # and it differs from the result with dhs_xmax = 1000, i.e. the cache was not used
    assert np.max(np.abs(direct.sed["thermal"] / first.sed["thermal"] - 1)) > 1e-8
