"""Temperature-dependent optical constants (DMS 2.1): file names, several components, the modes
'nearest' and 'linear', and the input of DMS 2.0.

The data sets are synthetic: the astrosilicate data with n and k changed linearly with the
temperature, written to a temporary optical-data directory.
"""

import os

import numpy as np
import pytest

from conftest import DATA_DIR, make_model, requires_dms

pytestmark = requires_dms

import dms  # noqa: E402

BASE_FILE = "astrosilicate_draine03_862"


def nk_at(temp):
    """Synthetic optical constants at temperature temp: n and k linear in T."""
    tab = np.loadtxt(os.path.join(DATA_DIR, BASE_FILE), comments="#")
    lam, n, k = tab[:, 0], tab[:, 1], tab[:, 2]
    return lam, n * (1.0 + temp / 1000.0), k * (0.5 + temp / 100.0)


def write_set(folder, name, temp):
    lam, n, k = nk_at(temp)
    os.makedirs(folder, exist_ok=True)
    np.savetxt(os.path.join(folder, name), np.column_stack([lam, n, k]), fmt="%.12e",
               header=f"synthetic test data, T = {temp} K")


def data_dir(tmp_path, sets, sub="tdep"):
    """Optical-data directory with astrosilicate and the given files {name: temperature}."""
    d = tmp_path / "data"
    os.makedirs(d / sub, exist_ok=True)
    src = os.path.join(DATA_DIR, BASE_FILE)
    if not os.path.exists(d / BASE_FILE):
        os.symlink(src, d / BASE_FILE)
    for name, temp in sets.items():
        write_set(str(d / sub), name, temp)
    return str(d) + "/"


def model(ddir, **dust):
    m = make_model(dust={"n_sizes": 10, **dust})
    m.update(run={"data_dir": ddir}, disk={"r_in": 20.0, "r_out": 120.0, "inclination": 60.0},
             wavelengths={"lambda": [1.6, 24.0, 70.0, 160.0]})
    return m


def close(a, b, rtol=1e-6):
    np.testing.assert_allclose(a.sed["thermal"], b.sed["thermal"], rtol=rtol)
    np.testing.assert_allclose(a.sed["scattered"], b.sed["scattered"], rtol=rtol)
    assert a.f_d == pytest.approx(b.f_d, rel=rtol)


def test_file_names(tmp_path):
    """<base>_<T>K[.ext]: integers, leading zeros, decimals, any or no extension; other names ignored."""
    sets = {"fo_10K.txt": 10.0, "fo_015K": 15.0, "fo_12.5K.dat": 12.5, "fo_200k.tab": 200.0}
    ddir = data_dir(tmp_path, sets)
    for decoy in ("fo_20K_old.txt", "fox_30K.txt", "fo_40.K", "fo_K.txt", "fo_50Kx", ".fo_60K.txt"):
        write_set(os.path.join(ddir, "tdep"), decoy, 99.0)
    found = dms.tdep_files("tdep/fo", ddir)
    assert [t for t, _ in found] == [10.0, 12.5, 15.0, 200.0]
    assert [f for _, f in found] == ["tdep/fo_10K.txt", "tdep/fo_12.5K.dat", "tdep/fo_015K",
                                     "tdep/fo_200k.tab"]
    res = model(ddir, optical_data="tdep/fo", tdep=True).run()
    assert list(res.species["T"]) == [10.0, 12.5, 15.0, 200.0]
    assert res.species["file"] == [f for _, f in found]
    assert res.species["is_data"].all()


def test_duplicate_temperature(tmp_path):
    """_010K and _10K are the same temperature: an error."""
    ddir = data_dir(tmp_path, {"fo_010K.txt": 10.0, "fo_10K.dat": 10.0, "fo_20K.txt": 20.0})
    with pytest.raises(ValueError, match="same temperature"):
        dms.tdep_files("tdep/fo", ddir)
    with pytest.raises(dms.DMSError, match="same temperature"):
        model(ddir, optical_data="tdep/fo", tdep=True).run()


def test_no_file_and_single_file(tmp_path):
    """No file: an error. One file: used without temperature dependence (same as a plain run)."""
    ddir = data_dir(tmp_path, {"one_50K.txt": 50.0})
    with pytest.raises(dms.DMSError, match="no files"):
        model(ddir, optical_data="tdep/none", tdep=True).run()
    a = model(ddir, optical_data="tdep/one", tdep=True).run()
    b = model(ddir, optical_data="tdep/one_50K.txt").run()
    close(a, b, rtol=1e-10)
    assert list(a.species["T"]) == [-1.0]


def test_both_inputs_rejected(tmp_path):
    ddir = data_dir(tmp_path, {"fo_10K.txt": 10.0, "fo_20K.txt": 20.0})
    m = model(ddir, optical_data="tdep/fo", tdep=True, tdep_component=1,
              tdep_temperatures=[10.0, 20.0], tdep_files=["tdep/fo_10K.txt", "tdep/fo_20K.txt"])
    with pytest.raises(dms.DMSError, match="not both"):
        m.run()


def test_input_of_version_2_0(tmp_path):
    """tdep_component/tdep_temperatures/tdep_files and tdep with the same files: same result."""
    sets = {"fo_30K.txt": 30.0, "fo_60K.txt": 60.0, "fo_100K.txt": 100.0}
    ddir = data_dir(tmp_path, sets)
    dust = {"n_components": 2, "bulk_density": [3.5, 3.0], "abundance": [1.0, 1.0]}
    old = model(ddir, optical_data=[BASE_FILE, "tdep/fo_30K.txt"], tdep_component=2,
                tdep_temperatures=[100.0, 30.0, 60.0],                   # any order
                tdep_files=["tdep/fo_100K.txt", "tdep/fo_30K.txt", "tdep/fo_60K.txt"], **dust).run()
    new = model(ddir, optical_data=[BASE_FILE, "tdep/fo"], tdep=[False, True], **dust).run()
    close(old, new, rtol=1e-10)
    np.testing.assert_array_equal(old.temperature["T"], new.temperature["T"])
    # the temperature dependence matters for this model (different from one fixed set)
    fixed = model(ddir, optical_data=[BASE_FILE, "tdep/fo_30K.txt"], **dust).run()
    assert np.max(np.abs(fixed.sed["thermal"] / new.sed["thermal"] - 1)) > 1e-3


def test_several_components_without_change(tmp_path):
    """Three temperature-dependent components whose sets are all identical: same result as the
    plain files, in both modes."""
    d = tmp_path / "data"
    os.makedirs(d / "tdep", exist_ok=True)
    os.symlink(os.path.join(DATA_DIR, BASE_FILE), d / BASE_FILE)
    for base in ("a", "b", "c"):
        for t in (20, 50, 90, 140):
            os.symlink(os.path.join(DATA_DIR, BASE_FILE), d / "tdep" / f"{base}_{t}K.txt")
    ddir = str(d) + "/"
    dust = {"n_components": 3, "bulk_density": [3.5, 3.0, 2.5], "abundance": [1.0, 2.0, 3.0]}
    plain = model(ddir, optical_data=[BASE_FILE] * 3, **dust).run()
    for mode in ("nearest", "linear"):
        res = model(ddir, optical_data=["tdep/a", "tdep/b", "tdep/c"], tdep=[True, True, True],
                    tdep_mode=mode, tdep_dt=10.0, **dust).run()
        close(plain, res, rtol=1e-10)


def test_linear_mode(tmp_path):
    """'linear': species between the sets with n, k interpolated in T, spaced by <= tdep_dt. With
    n, k linear in T, interpolating 20/60/100 K equals data files at 40 and 80 K."""
    coarse = data_dir(tmp_path, {"fo_20K.txt": 20.0, "fo_60K.txt": 60.0, "fo_100K.txt": 100.0})
    fine = data_dir(tmp_path / "f", {f"fo_{t}K.txt": float(t) for t in (20, 40, 60, 80, 100)})
    lin = model(coarse, optical_data="tdep/fo", tdep=True, tdep_mode="linear", tdep_dt=20.0).run()
    np.testing.assert_allclose(lin.species["T"], [20, 40, 60, 80, 100])
    assert list(lin.species["is_data"]) == [True, False, True, False, True]
    near = model(fine, optical_data="tdep/fo", tdep=True, tdep_mode="nearest").run()
    close(lin, near, rtol=1e-8)
    # step: 40 K intervals with tdep_dt = 15 K -> 3 steps of 13.3 K each
    lin15 = model(coarse, optical_data="tdep/fo", tdep=True, tdep_mode="linear", tdep_dt=15.0).run()
    assert np.max(np.diff(lin15.species["T"])) <= 15.0 + 1e-9
    assert len(lin15.species["T"]) == 7
    with pytest.raises(dms.DMSError, match="tdep_mode"):
        model(coarse, optical_data="tdep/fo", tdep=True, tdep_mode="cubic").run()


def test_temperature_outside_the_sets(tmp_path):
    """Grains colder than the coldest set: the coldest set is used and a warning is written."""
    ddir = data_dir(tmp_path, {"fo_300K.txt": 300.0, "fo_400K.txt": 400.0})
    m = model(ddir, optical_data="tdep/fo", tdep=True)
    m.update(run={"output_dir": str(tmp_path / "out") + "/", "output_name": "t", "write_files": True})
    a = m.run(write_files=True)
    b = model(ddir, optical_data="tdep/fo_300K.txt").run()
    close(a, b, rtol=1e-10)
    summary = open(tmp_path / "out" / "summary.t.txt").read()
    assert "below the coldest set" in summary
    assert "300.00 K  tdep/fo_300K.txt" in summary


def test_blowout_and_sublimation_use_the_local_set(tmp_path):
    """beta of a grain comes from the set used at its distance: in a disk that spans both sets the
    reported (largest) beta equals the larger of the two single-set values. Sublimation of a
    temperature-dependent component forms its residue only where the grains are hot enough."""
    ddir = data_dir(tmp_path, {"fo_40K.txt": 40.0, "fo_400K.txt": 400.0})
    geo = {"r_in": 3.0, "r_out": 150.0}
    star = {"teff": 8590.0, "radius": 1.84, "mass": 1.92}
    def run(**dust):
        m = model(ddir, a_min=0.3, **dust)
        m.update(disk=geo, star=star)
        return m.run()
    both = run(optical_data="tdep/fo", tdep=True)
    used = set(np.round(both.species["T"]).astype(int))
    assert used == {40, 400}
    cold = run(optical_data="tdep/fo_40K.txt")
    hot = run(optical_data="tdep/fo_400K.txt")
    np.testing.assert_allclose(both.grains["beta"][:, 0],
                               np.maximum(cold.grains["beta"][:, 0], hot.grains["beta"][:, 0]), rtol=1e-12)
    T = both.temperature["T"][:, 0, :]               # [n_dist, n_size]
    assert T.max() > 220.0 and T.min() < 220.0      # both sets are used somewhere
    sub = run(optical_data=["tdep/fo", BASE_FILE], n_components=2, bulk_density=[3.0, 3.5],
              abundance=[1.0, 0.0], tdep=[True, False], sublimation=True,
              t_sublimation=[150.0, -1.0], sublimation_residue=[2, 0])
    assert sub.grains["number"][:, 1].sum() > 0                 # residue formed inside the radius
    none = run(optical_data=["tdep/fo", BASE_FILE], n_components=2, bulk_density=[3.0, 3.5],
               abundance=[1.0, 0.0], tdep=[True, False], sublimation=True,
               t_sublimation=[2000.0, -1.0], sublimation_residue=[2, 0])
    assert none.grains["number"][:, 1].sum() == 0
