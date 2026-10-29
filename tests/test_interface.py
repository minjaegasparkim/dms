"""Input handling, error handling, Python interface and command-line program."""

import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

import dms
from conftest import DATA_DIR, DMS_EXE, make_model, requires_dms, requires_exe


@requires_dms
def test_invalid_input_raises_and_session_continues():
    with pytest.raises(dms.DMSError):
        make_model(disk={"r_in": 200.0, "r_out": 100.0}).run()
    with pytest.raises(dms.DMSError):
        make_model(dust={"optical_data": "no_such_file"}).run()
    res = make_model().run()                           # the next run works normally
    assert res.f_d > 0


@requires_dms
def test_repeated_runs_identical():
    r1 = make_model().run()
    r2 = make_model().run()
    # identical up to the summation order of the OpenMP threads
    assert r1.f_d == pytest.approx(r2.f_d, rel=1e-12)
    np.testing.assert_allclose(r1.sed["total"], r2.sed["total"], rtol=1e-10)


@requires_dms
def test_namelist_roundtrip(tmp_path):
    m = make_model(disk={"eccentricity": 0.2, "omega": 30.0}, dust={"blowout_release": "pericentre"})
    f = tmp_path / "model.nml"
    m.write(str(f))
    m2 = dms.Model.from_file(str(f))
    assert m2.run().f_d == pytest.approx(m.run().f_d, rel=1e-12)


@requires_dms
def test_result_contents():
    res = make_model(run={"mode": "both"}).run()
    nl = len(res.wavelength)
    assert res.image("total").shape == (nl, 41, 41)
    assert set(res.sed) >= {"thermal", "scattered", "Q", "U", "polarised", "star", "total"}
    np.testing.assert_allclose(res.sed["total"], res.sed["thermal"] + res.sed["scattered"]
                               + res.sed["star"], rtol=1e-10)
    assert res.tau_radial < 0.1


@requires_exe
@requires_dms
def test_command_line_matches_python(tmp_path):
    m = make_model(run={"mode": "sed", "output_dir": str(tmp_path / "out"), "output_name": "t",
                        "write_files": True})
    nml = tmp_path / "model.nml"
    m.write(str(nml))
    p = subprocess.run([DMS_EXE, str(nml)], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    sed = np.loadtxt(tmp_path / "out" / "SED.t.dat")
    res = m.run()
    np.testing.assert_allclose(sed[:, 1], res.sed["thermal"], rtol=1e-5)


@requires_exe
def test_command_line_error_exit_code(tmp_path):
    nml = tmp_path / "bad.nml"
    nml.write_text("&disk\n  r_in = 100., r_out = 50.\n/\n")
    p = subprocess.run([DMS_EXE, str(nml)], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode != 0


@requires_exe
def test_help():
    p = subprocess.run([DMS_EXE, "--help"], capture_output=True, text=True)
    assert p.returncode == 0 and "usage" in p.stdout


def test_data_directory_present():
    assert os.path.isfile(os.path.join(DATA_DIR, "astrosilicate_draine03_862"))


@requires_dms
def test_quickstart_example(tmp_path):
    """examples/quickstart.py (the example on the front page) runs and writes its figure."""
    pytest.importorskip("matplotlib")
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples",
                          "quickstart.py")
    env = dict(os.environ, MPLBACKEND="Agg")
    p = subprocess.run([sys.executable, script], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "f_d = " in p.stdout
    assert (tmp_path / "quickstart.png").is_file()


def test_size_grid_warning():
    """Fewer than 25 size bins per decade give a SizeGridWarning; the default (100 bins for
    1-1000 micron) does not."""
    with pytest.warns(dms.SizeGridWarning):
        dms.model.check_size_grid({"n_sizes": 30})
    with warnings.catch_warnings():
        warnings.simplefilter("error", dms.SizeGridWarning)
        dms.model.check_size_grid({})
        dms.model.check_size_grid({"size_distribution": "single", "n_sizes": 1})


@requires_dms
def test_default_number_of_size_bins():
    m = make_model()
    del m.params["dust"]["n_sizes"]
    res = m.run()
    assert len(res.grains["a"]) == 100


@requires_dms
@pytest.mark.skipif(not os.environ.get("DMS_TEST_SLOW"), reason="about two minutes; set DMS_TEST_SLOW=1")
def test_two_belts_observed_example(tmp_path):
    """examples/two_belts_observed.py (Example 2 of the README) runs and writes its figure."""
    pytest.importorskip("matplotlib")
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples",
                          "two_belts_observed.py")
    env = dict(os.environ, MPLBACKEND="Agg")
    p = subprocess.run([sys.executable, script], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "blow-out threshold       beta >= 0.44 ... 0.50" in p.stdout
    assert (tmp_path / "two_belts_observed.png").is_file()


@requires_exe
@requires_dms
def test_parameter_study_example(tmp_path):
    """examples/parameter_study.py (Example 3 of the README): 12 models with dms.batch, then a
    second call that finds all results and computes nothing."""
    pytest.importorskip("matplotlib")
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples",
                          "parameter_study.py")
    env = dict(os.environ, MPLBACKEND="Agg")
    for expected in ("12 models, 12 to compute", "12 models, 0 to compute"):
        p = subprocess.run([sys.executable, script], cwd=tmp_path, capture_output=True, text=True, env=env)
        assert p.returncode == 0, p.stdout + p.stderr
        assert expected in p.stdout
    assert "ice50_G_q3.5" in p.stdout
    assert (tmp_path / "parameter_study.png").is_file()
