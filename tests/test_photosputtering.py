"""Ice loss by photosputtering (DMS 2.2): the switch, the ice-loss radius against the analytic
expression, the reference times, the residue, the UV spectrum file and the input checks.

Model: an A star (8590 K, 1.84 R_sun, 1.92 M_sun) and a disk from 5 to 300 au of an ice-silicate
mixture (component 1, ice mass fraction 0.3) that becomes astronomical silicate (component 2)
when its ice is gone.
"""

import os

import numpy as np
import pytest

from conftest import make_model, requires_dms

pytestmark = requires_dms

import dms  # noqa: E402

trapezoid = getattr(np, "trapezoid", None) or np.trapz   # np.trapezoid from NumPy 2.0 on

M_H2O = 18.01528 * 1.66053906660e-27          # [kg]
YEAR = 365.25 * 86400.0
AU = 149597870700.0
STAR = {"teff": 8590.0, "radius": 1.84, "mass": 1.92, "distance": 10.0}
RHO_W = 0.3 * 2300.0                           # ice mass per grain volume [kg/m^3]


def model(**dust):
    d = {"n_components": 2, "optical_data": ["ci50_sil50_BR", "astrosilicate_draine03_862"],
         "bulk_density": [2.3, 3.5], "abundance": [1.0, 0.0], "a_min": 1.0, "n_sizes": 12,
         "photosputtering": True, "ps_yield": [1e-3, 0.0], "ps_ice_fraction": [0.3, 1.0],
         "ps_residue": [2, 0]}
    d.update(dust)
    return make_model(star=STAR, disk={"r_in": 5.0, "r_out": 300.0}, dust=d)


def test_switch_off_changes_nothing():
    """photosputtering = .false. with the other inputs set gives the result of a run without them."""
    plain = make_model(star=STAR, disk={"r_in": 5.0, "r_out": 300.0},
                       dust={"n_components": 2, "optical_data": ["ci50_sil50_BR", "astrosilicate_draine03_862"],
                             "bulk_density": [2.3, 3.5], "abundance": [1.0, 0.0], "a_min": 1.0,
                             "n_sizes": 12}).run()
    off = model(photosputtering=False, ps_time="age", ps_age=1e3).run()
    np.testing.assert_allclose(off.sed["thermal"], plain.sed["thermal"], rtol=1e-10)
    np.testing.assert_allclose(off.sed["scattered"], plain.sed["scattered"], rtol=1e-10)
    assert off.photosputtering is None


def test_age_mode_against_analytic_radius():
    """Q = 1 and t_ref = age: the ice is lost inside r = sqrt(N_1au m Y t / (4 rho_w a)) au."""
    res = model(ps_time="age", ps_age=1e3, ps_absorption="unity").run()
    ps = res.photosputtering
    a = res.grains["a"]
    n1 = ps["photons_1au"] * 1e4
    r_an = np.sqrt(n1 * M_H2O * 1e-3 * 1e3 * YEAR / (4.0 * RHO_W * a * 1e-6))
    inside = (r_an > 10.0) & (r_an < 250.0)
    assert inside.sum() >= 5
    np.testing.assert_allclose(ps["r_ice_loss"][inside, 0], r_an[inside], rtol=1e-3)
    assert np.all(ps["r_ice_loss"][:, 1] == -1.0)              # the residue has no photosputtering
    # photon flux of the black body in 5.1 - 13.6 eV at 1 au
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    lam = np.geomspace(0.0912e-6, 0.2431e-6, 20001)
    B = 2 * h * c**2 / lam**5 / np.expm1(h * c / (lam * k * STAR["teff"]))
    ref = trapezoid(np.pi * B * (STAR["radius"] * 6.957e8 / AU)**2 * lam / (h * c), lam) * 1e-4
    assert ps["photons_1au"] == pytest.approx(ref, rel=1e-4)


def test_residue_and_vanishing_grains():
    """The grains that lose their ice become residue grains (ps_residue = 2) or vanish (0)."""
    kept = model(ps_time="age", ps_age=1e3).run()
    gone = model(ps_time="age", ps_age=1e3, ps_residue=[0, 0]).run()
    assert kept.grains["number"][:, 1].sum() > 0
    assert gone.grains["number"][:, 1].sum() == 0
    # grains that are unbound as icy grains are blown out as they are: no residue of those sizes
    unbound = kept.grains["beta"][:, 0] >= 0.5
    assert unbound.any()
    assert np.all(kept.grains["number"][unbound, 1] == 0)
    np.testing.assert_array_equal(kept.photosputtering["r_ice_loss"], gone.photosputtering["r_ice_loss"])


def test_pr_mode_does_not_depend_on_distance():
    """t_ps and t_PR both grow as r^2: a grain size loses its ice everywhere or nowhere, decided by
    4 rho_w a 2 beta G M / (m Y N_1au au^2 c) < 1."""
    res = model(ps_time="pr", ps_absorption="unity").run()
    rl = res.photosputtering["r_ice_loss"][:, 0]
    rmax = rl.max()
    assert np.all((rl == 0.0) | (rl == rmax))
    beta = res.grains["beta"][:, 0]
    a = res.grains["a"] * 1e-6
    n1 = res.photosputtering["photons_1au"] * 1e4
    x = 4 * RHO_W * a * 2 * beta * 6.67430e-11 * STAR["mass"] * 1.98840987e30 / (M_H2O * 1e-3 * n1 * AU**2 * 2.99792458e8)
    clear = (beta < 0.5) & (np.abs(np.log(x)) > 0.05)       # unbound grains are not tested
    assert clear.sum() >= 8
    assert np.all(rl[beta >= 0.5] == 0.0)
    np.testing.assert_array_equal((rl > 0)[clear], (x < 1)[clear])


def test_collision_mode_and_shortest():
    """A larger optical depth (shorter collisional lifetime) keeps more ice; 'shortest' takes the
    smallest of the reference times, so its ice-loss radius is the smallest of the single ones."""
    with pytest.raises(dms.DMSError, match="ps_tau"):
        model(ps_time="collision").run()
    lo = model(ps_time="collision", ps_tau=1e-5).run().photosputtering["r_ice_loss"][:, 0]
    hi = model(ps_time="collision", ps_tau=1e-2).run().photosputtering["r_ice_loss"][:, 0]
    assert np.all(hi <= lo) and np.any(hi < lo)
    age = model(ps_time="age", ps_age=3e3).run().photosputtering["r_ice_loss"][:, 0]
    pr = model(ps_time="pr").run().photosputtering["r_ice_loss"][:, 0]
    short = model(ps_time="shortest", ps_tau=1e-5, ps_age=3e3).run().photosputtering["r_ice_loss"][:, 0]
    np.testing.assert_array_equal(short, np.minimum(np.minimum(lo, age), pr))


def test_uv_spectrum_file(tmp_path):
    """A UV spectrum with ten times the black-body flux: ten times the photons, sqrt(10) times the
    ice-loss radius (age mode)."""
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    lam = np.geomspace(0.08, 0.30, 400)                                   # micron
    B = 2 * h * c**2 / (lam * 1e-6)**5 / np.expm1(h * c / (lam * 1e-6 * k * STAR["teff"]))
    d = STAR["distance"] * 648000.0 / np.pi * AU
    fnu = 10 * np.pi * B * (STAR["radius"] * 6.957e8 / d)**2 * (lam * 1e-6)**2 / c * 1e26   # Jy
    f = tmp_path / "uv.txt"
    np.savetxt(f, np.column_stack([lam, fnu]), header="lambda [micron]  F_nu [Jy]")
    base = model(ps_time="age", ps_age=100.0, ps_absorption="unity").run().photosputtering
    uv = model(ps_time="age", ps_age=100.0, ps_absorption="unity", ps_uv_file=str(f)).run().photosputtering
    assert uv["photons_1au"] == pytest.approx(10 * base["photons_1au"], rel=2e-3)
    ok = (base["r_ice_loss"][:, 0] > 10) & (uv["r_ice_loss"][:, 0] < 250)
    assert ok.sum() >= 3
    np.testing.assert_allclose(uv["r_ice_loss"][ok, 0], np.sqrt(10) * base["r_ice_loss"][ok, 0], rtol=2e-3)


def test_interstellar_field_and_summary(tmp_path):
    """Interstellar UV adds ice loss in the outer disk; the summary lists the settings."""
    m = model(ps_time="age", ps_age=1e3, ps_isrf=1e10)
    m.update(run={"output_dir": str(tmp_path) + "/", "output_name": "ps", "write_files": True})
    res = m.run(write_files=True)
    base = model(ps_time="age", ps_age=1e3).run()
    assert np.all(res.photosputtering["r_ice_loss"][:, 0] >= base.photosputtering["r_ice_loss"][:, 0])
    assert np.any(res.photosputtering["r_ice_loss"][:, 0] > base.photosputtering["r_ice_loss"][:, 0])
    text = open(os.path.join(tmp_path, "summary.ps.txt")).read()
    assert "photosputtering: ice lost where the sputtering time < age" in text
    assert "interstellar UV" in text and "ice lost out to" in text


@pytest.mark.parametrize("dust, msg", [
    ({"ps_yield": [0.0, 0.0]}, "ps_yield"),
    ({"ps_yield": [-1e-3, 0.0]}, "ps_yield"),
    ({"ps_ice_fraction": [0.0, 1.0]}, "ps_ice_fraction"),
    ({"ps_residue": [1, 0]}, "ps_residue"),
    ({"ps_yield": [1e-3, 1e-3]}, "ps_residue"),
    ({"ps_time": "forever"}, "ps_time"),
    ({"ps_time": "age"}, "ps_age"),
    ({"ps_absorption": "half"}, "ps_absorption"),
    ({"ps_band": [0.3, 0.1]}, "ps_band"),
])
def test_invalid_input(dust, msg):
    with pytest.raises(dms.DMSError, match=msg):
        model(**dust).run()


# ------------------------------------------------------------------------------------------------
# ps_orbit = 'release': UV averaged over the orbit of the grain after its release
# ------------------------------------------------------------------------------------------------
def orbit_factor(beta, e=0.0, release="circular"):
    """r^2 <r^-2> on the orbit after release at r from a parent orbit (a_p, e) (in units of r)."""
    a_p = {"circular": 1.0, "pericentre": 1.0 / (1.0 - e), "apocentre": 1.0 / (1.0 + e)}[release]
    if release == "circular":
        e = 0.0
    a_g = (1.0 - beta) / (1.0 / a_p - 2.0 * beta)
    ome2 = a_p * (1.0 - e**2) / ((1.0 - beta) * a_g)
    return 1.0 / (a_g**2 * np.sqrt(ome2)), a_g, ome2


def test_release_orbit_circular_parents():
    """Circular parents: N_abs is scaled by (1-2 beta)^(3/2) / (1-beta), so in age mode the
    ice-loss radius by its square root; ps_orbit = 'local' is the default."""
    kw = dict(ps_time="age", ps_age=1e3, ps_absorption="unity")
    default = model(**kw).run()
    local = model(ps_orbit="local", **kw).run()
    np.testing.assert_array_equal(default.photosputtering["r_ice_loss"], local.photosputtering["r_ice_loss"])
    rel = model(ps_orbit="release", **kw).run()
    beta = rel.grains["beta"][:, 0]
    bnd = np.minimum(beta, 0.49)                                          # unbound grains: blown out
    f = (1 - 2 * bnd)**1.5 / (1 - bnd)
    np.testing.assert_allclose(f, orbit_factor(bnd)[0], rtol=1e-12)
    rl, rr = local.photosputtering["r_ice_loss"][:, 0], rel.photosputtering["r_ice_loss"][:, 0]
    ok = (beta < 0.5) & (rl > 10) & (rl < 250) & (rr > 10)
    assert ok.sum() >= 5 and beta[ok].max() > 0.05
    np.testing.assert_allclose(rr[ok], rl[ok] * np.sqrt(f[ok]), rtol=1e-3)


def eccentric(release="pericentre", e=0.3, bins=32, mode="sed", belts=None, **dust):
    """Eccentric belt (a = 5-300 au) of the icy mixture, photosputtering on release orbits."""
    d = {"n_components": 2, "optical_data": ["ci50_sil50_BR", "astrosilicate_draine03_862"],
         "bulk_density": [2.3, 3.5], "abundance": [1.0, 0.0], "a_min": 1.0, "n_sizes": 12,
         "blowout_release": release, "photosputtering": True, "ps_yield": [1e-3, 0.0],
         "ps_ice_fraction": [0.3, 1.0], "ps_residue": [2, 0], "ps_orbit": "release",
         "ps_orbit_bins": bins}
    d.update(dust)
    disk = {"r_in": 5.0, "r_out": 300.0, "eccentricity": e, "omega": 30.0, "inclination": 50.0}
    if belts:
        disk = belts
    m = make_model(star=STAR, disk=disk, dust=d)
    m.update(run={"mode": mode})
    return m


@pytest.mark.parametrize("release", ["pericentre", "apocentre"])
def test_eccentric_belt_ice_state_per_orbit(release):
    """Eccentric belt (e = 0.3): the ice state is decided per orbit. In age mode with Q = 1 an orbit
    loses the ice of a grain if 4 rho_w a_g^2 sqrt(1 - e_g^2) / (m Y N_1au au^2) < age for the orbit
    after release at pericentre or apocentre of the parent orbit (a, e): inside a_crit ∝ 1/s.
    a_ice_loss is the centre of the last orbit bin inside a_crit (200 bins: 2 % wide)."""
    res = eccentric(release, bins=200, ps_time="age", ps_age=1e3, ps_absorption="unity").run()
    ps = res.photosputtering
    assert ps["orbit_bins"] == 200
    assert np.all(ps["r_ice_loss"][:, 0] == 0.0)          # no circular belt
    beta = res.grains["beta"][:, 0]
    s = res.grains["a"] * 1e-6
    n1 = ps["photons_1au"] * 1e4
    e = 0.3
    x = 1 - e if release == "pericentre" else 1 + e      # release distance / a
    bound = beta < x / 2
    beta = np.where(bound, beta, 0.0)                     # unbound sizes are not tested
    a_g = (1 - beta) / (1 - 2 * beta / x)                 # in units of a
    ome2 = (1 - e**2) / ((1 - beta) * a_g)
    # t_ps = 4 rho_w s a^2 a_g^2 sqrt(ome2) / (m Y N_1au au^2) < t  ->  a < a_crit
    a_crit = np.sqrt(1e3 * YEAR * M_H2O * 1e-3 * n1 / (4 * RHO_W * s * a_g**2 * np.sqrt(ome2))) * 1.0   # [au]
    width = (300.0 / 5.0)**(1 / 200)
    ok = bound & (a_crit > 6) & (a_crit < 280)
    assert ok.sum() >= 4
    al = ps["a_ice_loss"][ok, 0]
    assert np.all(al <= a_crit[ok]) and np.all(al * width > a_crit[ok] / width)


def test_eccentric_belt_without_ice_loss_matches_standard_run():
    """With a negligible yield no grain loses its ice: the orbit cells (two belts, one of them
    eccentric) must reproduce the run without photosputtering in SED, images, f_d, optical depth
    and grain numbers."""
    belts = {"n_belts": 2, "r_in": [5.0, 80.0], "r_out": [40.0, 150.0], "eccentricity": [0.3, 0.0],
             "omega": 30.0, "belt_mass_fraction": [0.4, 0.6], "inclination": 50.0}
    kw = dict(belts=belts, sublimation=True, t_sublimation=[110.0, -1.0], sublimation_residue=[2, 0])
    on = eccentric(mode="both", ps_yield=[1e-30, 0.0], ps_time="pr", **kw).run()
    off = eccentric(mode="both", photosputtering=False, **kw).run()
    assert on.photosputtering["orbit_bins"] == 32
    for k in ("thermal", "scattered"):
        np.testing.assert_allclose(on.sed[k], off.sed[k], rtol=1e-10)
    for k in off.images:
        a, b = off._images[k], on._images[k]
        assert np.nanmax(np.abs(b - a)) <= 1e-10 * np.nanmax(np.abs(a))
    assert on.f_d == pytest.approx(off.f_d, rel=1e-10)
    assert on.tau_peak == pytest.approx(off.tau_peak, rel=1e-10)
    np.testing.assert_allclose(on.grains["number"], off.grains["number"], rtol=1e-10, atol=1e-30)


def test_eccentric_belt_small_e_matches_circular_belt():
    """e -> 0: the orbit bins reproduce the circular belt, where the ice is decided per distance
    (SED within the 2 % width of the 200 orbit bins)."""
    kw = dict(ps_time="age", ps_age=1e3)
    ecc = eccentric("pericentre", e=1e-4, bins=200, **kw).run()
    circ = eccentric("pericentre", e=0.0, **kw).run()
    assert "a_ice_loss" in ecc.photosputtering and "a_ice_loss" not in circ.photosputtering
    for k in ("thermal", "scattered"):
        np.testing.assert_allclose(ecc.sed[k], circ.sed[k], rtol=2e-2)
    lost = circ.photosputtering["r_ice_loss"][:, 0]
    a_l = ecc.photosputtering["a_ice_loss"][:, 0]
    inside = (lost > 6) & (lost < 280)
    assert inside.sum() >= 2
    np.testing.assert_allclose(a_l[inside], lost[inside], rtol=0.05)


def test_eccentric_belt_ice_loss_differs_from_local():
    """The point of the orbit bins: in an eccentric belt the local test (per distance, as for a
    circular belt) and the orbit test give different SEDs."""
    kw = dict(ps_time="age", ps_age=1e3, ps_absorption="unity")
    orb = eccentric("pericentre", **kw).run()
    loc = eccentric("pericentre", ps_orbit="local", **kw).run()
    assert "a_ice_loss" in orb.photosputtering and "a_ice_loss" not in loc.photosputtering
    assert np.max(np.abs(orb.sed["thermal"] / loc.sed["thermal"] - 1)) > 1e-3


def test_release_orbit_pr_time():
    """P-R mode with release orbits: t_ps / f against a^2 c (1-e^2)^(3/2) / (beta G M (2 + 3 e^2));
    both grow as r^2, so each size again loses its ice everywhere or nowhere."""
    res = model(ps_time="pr", ps_absorption="unity", ps_orbit="release").run()
    rl = res.photosputtering["r_ice_loss"][:, 0]
    beta = res.grains["beta"][:, 0]
    bound = beta < 0.5
    assert np.all((rl[bound] == 0.0) | (rl[bound] == rl.max()))
    f, a_g, ome2 = orbit_factor(np.where(bound, beta, 0.1))
    a = res.grains["a"] * 1e-6
    n1 = res.photosputtering["photons_1au"] * 1e4
    GM = 6.67430e-11 * STAR["mass"] * 1.98840987e30
    # per r^2: t_ps r^-2 = 4 rho_w a / (m Y N_1au au^2 f); t_PR r^-2 = a_g^2 c ome2^1.5 / (beta GM (2 + 3 e^2))
    tps = 4 * RHO_W * a / (M_H2O * 1e-3 * n1 * AU**2 * f)
    tpr = a_g**2 * 2.99792458e8 * ome2**1.5 / (beta * GM * (2 + 3 * (1 - ome2)))
    x = tps / tpr
    clear = bound & (np.abs(np.log(x)) > 0.05)
    assert clear.sum() >= 6
    np.testing.assert_array_equal((rl > 0)[clear], (x < 1)[clear])


def test_invalid_orbit():
    with pytest.raises(dms.DMSError, match="ps_orbit"):
        model(ps_orbit="average").run()
    with pytest.raises(dms.DMSError, match="ps_orbit_bins"):
        model(ps_orbit="release", ps_orbit_bins=0).run()


def test_eccentric_belt_energy_balance_and_images():
    """With ice lost on part of the orbits: the integrated thermal emission equals f_d L*, and the
    SED summed from the images agrees with the direct quadrature as well as without photosputtering."""
    from test_physics import thermal_luminosity_ratio
    kw = dict(ps_time="age", ps_age=1e3, ps_absorption="unity")
    m = eccentric("pericentre", **kw)
    m.update(run={"emission": "thermal"})
    m.params["wavelengths"] = {"lambda_min": 0.5, "lambda_max": 2000.0, "n_lambda": 120}
    res = m.run()
    assert 0 < np.count_nonzero(res.photosputtering["a_ice_loss"][:, 0] > 6) < 12
    assert thermal_luminosity_ratio(res, 10.0) / res.f_d == pytest.approx(1.0, abs=5e-3)
    # images against the direct quadrature: the same agreement as without photosputtering (the
    # wide eccentric disk on 61 pixels limits it to about 1 % either way)
    ratio = {}
    for ps in (True, False):
        sed = eccentric("pericentre", photosputtering=ps, **kw)
        sed.update(run={"emission": "thermal"})
        img = eccentric("pericentre", mode="both", photosputtering=ps, **kw)
        img.update(run={"emission": "thermal"}, image={"n_pixels": 61, "n_split_max": 3})
        s, i = sed.run(), img.run()
        ratio[ps] = (i.sed["thermal"] / s.sed["thermal"], i.f_d / s.f_d)
    np.testing.assert_allclose(ratio[True][0], ratio[False][0], atol=1e-3)
    assert ratio[True][1] == pytest.approx(ratio[False][1], abs=1e-3)
