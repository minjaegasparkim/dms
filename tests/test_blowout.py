"""Blow-out threshold (&dust blowout, blowout_beta, blowout_release)."""

import numpy as np
import pytest

import dms
from conftest import make_model, requires_dms

pytestmark = requires_dms

SMALL = {"a_min": 0.1, "n_sizes": 40}


def surviving_beta(res):
    beta, number = res.grains["beta"][:, 0], res.grains["number"][:, 0]
    return beta[number > 0], beta[number == 0]


def test_default_is_circular():
    """Without new inputs the threshold is 0.5, as in earlier versions."""
    r0 = make_model(dust=SMALL).run()
    r1 = make_model(dust=dict(SMALL, blowout_release="circular")).run()
    assert r0.blowout_beta == (0.5, 0.5)
    assert r1.f_d == pytest.approx(r0.f_d, rel=1e-12)
    kept, removed = surviving_beta(r0)
    assert kept.max() < 0.5 <= removed.min()


@pytest.mark.parametrize("release, crit", [("pericentre", 0.35), ("apocentre", 0.65)])
def test_release_in_eccentric_belt(release, crit):
    """e = 0.3: threshold (1 - e)/2 at pericentre, (1 + e)/2 at apocentre."""
    res = make_model(disk={"eccentricity": 0.3}, dust=dict(SMALL, blowout_release=release)).run()
    assert res.blowout_beta == pytest.approx((crit, crit))
    kept, removed = surviving_beta(res)
    assert kept.max() < crit <= removed.min()


def test_pericentre_release_removes_more_grains():
    circ = make_model(disk={"eccentricity": 0.3}, dust=SMALL).run()
    peri = make_model(disk={"eccentricity": 0.3}, dust=dict(SMALL, blowout_release="pericentre")).run()
    assert peri.grains["number"][:, 0].astype(bool).sum() < circ.grains["number"][:, 0].astype(bool).sum()


def test_circular_belt_pericentre_equals_circular():
    """For e = 0 all release points give 1/2."""
    r0 = make_model(dust=SMALL).run()
    r1 = make_model(dust=dict(SMALL, blowout_release="pericentre")).run()
    assert r1.blowout_beta == (0.5, 0.5)
    assert r1.f_d == pytest.approx(r0.f_d, rel=1e-12)


def test_user_threshold():
    res = make_model(dust=dict(SMALL, blowout_beta=0.2)).run()
    assert res.blowout_beta == pytest.approx((0.2, 0.2))
    kept, removed = surviving_beta(res)
    assert kept.max() < 0.2 <= removed.min()


def test_threshold_per_belt():
    """Two belts with e = 0 and 0.4, release at pericentre: thresholds 0.5 and 0.3."""
    res = make_model(disk={"n_belts": 2, "density_law": ["powerlaw", "powerlaw"], "alpha": [1.0, 1.0],
                           "r_in": [10.0, 80.0], "r_out": [20.0, 100.0], "eccentricity": [0.0, 0.4],
                           "omega": [0.0, 0.0], "belt_mass_fraction": [0.5, 0.5]},
                     dust=dict(SMALL, blowout_release="pericentre")).run()
    assert res.blowout_beta == pytest.approx((0.3, 0.5))


def test_blowout_off():
    res = make_model(dust=dict(SMALL, blowout=False)).run()
    assert res.blowout_beta == (0.0, 0.0)
    assert np.all(res.grains["number"][:, 0] > 0)


def test_conflicting_inputs_raise():
    with pytest.raises(dms.DMSError, match="blowout"):
        make_model(dust=dict(SMALL, blowout_beta=0.3, blowout_release="pericentre")).run()
    with pytest.raises(dms.DMSError, match="blowout_release"):
        make_model(dust=dict(SMALL, blowout_release="periastron")).run()
