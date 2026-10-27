"""Physical consistency checks of DMS (independent of stored reference values)."""

import numpy as np
import pytest

from conftest import make_model, requires_dms

pytestmark = requires_dms

C = 2.99792458e8
AU = 149597870700.0
PC = 648000.0 / np.pi * AU
LSUN = 3.828e26


def thermal_luminosity_ratio(res, distance_pc):
    """Integrated thermal SED divided by the stellar luminosity."""
    nu = C / (res.wavelength * 1e-6)
    order = np.argsort(nu)
    l_th = 4 * np.pi * (distance_pc * PC) ** 2 * np.trapezoid(res.sed["thermal"][order] * 1e-26, nu[order])
    return l_th / (res.luminosity_lsun * LSUN)


def test_energy_balance_two_eccentric_belts():
    """The integrated thermal emission equals f_d L* (energy conservation)."""
    m = make_model(
        run={"emission": "thermal"},
        disk={"n_belts": 2, "density_law": ["powerlaw", "gaussian"], "alpha": [0.0, 0.0],
              "mu": [0.0, 140.0], "sigma": [0.0, 8.0], "r_in": [5.0, 110.0], "r_out": [12.0, 170.0],
              "eccentricity": [0.0, 0.12], "omega": [0.0, 40.0], "belt_mass_fraction": [0.05, 0.95],
              "vertical": "gaussian", "h_over_r": 0.03},
    )
    m.params["wavelengths"] = {"lambda_min": 0.5, "lambda_max": 2000.0, "n_lambda": 120}
    res = m.run()
    assert thermal_luminosity_ratio(res, 10.0) / res.f_d == pytest.approx(1.0, abs=5e-3)


def test_image_and_quadrature_sed_agree():
    """Thermal SED from the summed images equals the direct quadrature."""
    sed = make_model(run={"emission": "thermal"}).run()
    img = make_model(run={"mode": "both", "emission": "thermal"},
                     image={"n_pixels": 61, "n_split_max": 3}).run()
    sig = sed.sed["thermal"] > 1e-6 * sed.sed["thermal"].max()     # skip the Wien tail at 1.6 um
    np.testing.assert_allclose(img.sed["thermal"][sig], sed.sed["thermal"][sig], rtol=1e-2)
    assert img.f_d == pytest.approx(sed.f_d, rel=1e-2)


def test_belt_mass_fractions():
    """Two belts with mass fractions 0.2 / 0.8: recovered from the vertical optical depth."""
    m = make_model(disk={"n_belts": 2, "density_law": ["powerlaw", "powerlaw"], "alpha": [1.0, 1.0],
                         "r_in": [10.0, 60.0], "r_out": [20.0, 100.0], "belt_mass_fraction": [0.2, 0.8]},
                   dust={"a_min": 5.0})                # no blow-out within the size range
    res = m.run()
    r, tau = res.tau_profile["r_au"], res.tau_profile["mean"]
    w = tau * 2 * np.pi * r * np.gradient(r)          # cross section per radial bin
    inner = w[r < 40].sum() / w.sum()
    assert inner == pytest.approx(0.2, abs=0.01)


def test_size_distribution_from_file(tmp_path):
    """A size distribution read from a file equals the same power law."""
    ref = make_model(run={"emission": "thermal"}).run()
    a, e = ref.grains["a"], ref.grains["edges"]
    q = 3.5
    w = (e[1:] ** (1 - q) - e[:-1] ** (1 - q)) / (1 - q)
    f = tmp_path / "sizes.txt"
    np.savetxt(f, np.c_[a, w])
    res = make_model(run={"emission": "thermal"},
                     dust={"size_distribution": "file", "size_file": str(f)}).run()
    assert res.f_d == pytest.approx(ref.f_d, rel=1e-6)
    np.testing.assert_allclose(res.sed["thermal"], ref.sed["thermal"], rtol=1e-6)


def test_polarisation_face_on_rayleigh():
    """Face-on ring of grains much smaller than the wavelength (Mie, x = 0.2): every grain is seen
    at a scattering angle near 90 deg, so the polarisation is close to 100 per cent, centro-symmetric
    (Q_phi = PI > 0) and U_phi = 0."""
    res = make_model(run={"mode": "image", "emission": "polarisation"},
                     disk={"inclination": 0.0}, wavelengths={"lambda": [1.6]},
                     dust={"size_distribution": "single", "a_min": 0.05, "blowout": False}).run()
    qphi, uphi, pi_ = (res.image(n)[0] for n in ("Qphi", "Uphi", "polarised"))
    # a pixel collects sub-volumes at slightly different position angles, so U_phi is not exactly
    # zero; the residual is of the order of (sub-volume / distance)^2
    assert qphi.sum() > 0
    assert np.abs(uphi).sum() < 2e-3 * qphi.sum()
    assert qphi.sum() == pytest.approx(pi_.sum(), rel=1e-3)
    deg = res.sed["polarised"][0] / res.sed["scattered"][0]
    assert deg == pytest.approx(0.0, abs=1e-6)         # face-on: Q and U cancel in the SED
    dop = pi_.sum() / res.image("scattered")[0].sum()
    assert dop > 0.97


def test_henyey_greenstein_polarisation_law():
    """HG phase function with p_max = 0.4, face-on ring with a 5 deg half opening angle: the
    scattering angles lie within 90 +- 5 deg, where p_max sin^2/(1 + cos^2) is between
    0.985 p_max and p_max."""
    res = make_model(run={"mode": "image", "emission": "polarisation"},
                     disk={"inclination": 0.0}, wavelengths={"lambda": [1.6]},
                     scattering={"phase_function": "hg", "g1": 0.5, "p_max": 0.4}).run()
    dop = res.image("polarised")[0].sum() / res.image("scattered")[0].sum()
    assert 0.985 * 0.4 < dop <= 0.4


def test_position_angle_rotates_image():
    """PA = 90 deg rotates the image of PA = 0 by 90 deg (counter-clockwise)."""
    kw = dict(run={"mode": "image", "emission": "thermal"}, wavelengths={"lambda": [70.0]})
    img0 = make_model(**kw).run().image("thermal")[0]
    kw["disk"] = {"position_angle": 90.0}
    img90 = make_model(**kw).run().image("thermal")[0]
    # arrays are [iy, ix] with iy increasing upwards; a counter-clockwise rotation on the sky
    rot = np.rot90(img0, k=1)
    np.testing.assert_allclose(img90, rot, rtol=1e-6, atol=1e-9 * img0.max())


def test_forward_scattering_on_near_side():
    """The near side (-X, facing the observer) is brighter in scattered light."""
    res = make_model(run={"mode": "image", "emission": "scattering"},
                     wavelengths={"lambda": [1.6]}, dust={"a_min": 2.0}).run()
    img = res.image("scattered")[0]
    n = img.shape[1]
    # PA = 0: major axis along y; near side at negative x (columns left of the centre)
    near = img[:, : n // 2].sum()
    far = img[:, n // 2 + 1:].sum()
    assert near > far


def test_eccentric_belt_extent():
    """Eccentric belt: dust between r_in (1 - e) and r_out (1 + e)."""
    e = 0.3
    res = make_model(disk={"eccentricity": e, "omega": 0.0, "r_in": 60.0, "r_out": 80.0}).run()
    r = res.tau_profile["r_au"]
    peri, apo = res.tau_profile["pericentre"], res.tau_profile["apocentre"]
    assert r[peri > 0].min() == pytest.approx(60.0 * (1 - e), rel=0.02)
    assert r[apo > 0].max() == pytest.approx(80.0 * (1 + e), rel=0.02)


def test_fractional_luminosity_scales_with_mass():
    """Optically thin: f_d and the SED are proportional to the dust mass."""
    r1 = make_model().run()
    r2 = make_model(disk={"dust_mass": 3.0e-8}).run()
    assert r2.f_d / r1.f_d == pytest.approx(3.0, rel=1e-10)
    np.testing.assert_allclose(r2.sed["thermal"], 3.0 * r1.sed["thermal"], rtol=1e-10)
