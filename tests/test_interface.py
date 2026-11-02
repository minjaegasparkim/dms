"""Input handling, error handling, Python interface and command-line program."""

import json
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

import dms
from conftest import DATA_DIR, DMS_EXE, EXAMPLES, REPO, make_model, requires_dms, requires_exe


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


# ------------------------------------------------------------------------------------------
# the example notebooks (examples/*.ipynb): their code cells are run as one script in a
# temporary directory, so that nothing is written into the repository
# ------------------------------------------------------------------------------------------

def run_notebook(name, tmp_path):
    pytest.importorskip("matplotlib")
    with open(os.path.join(EXAMPLES, name)) as f:
        nb = json.load(f)
    code = "\n\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    script = tmp_path / (name[:-6] + ".py")
    script.write_text(code)
    env = dict(os.environ, MPLBACKEND="Agg", DMS_REPO=REPO)
    p = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True, text=True, env=env)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout


@requires_dms
def test_notebook_0_getting_started(tmp_path):
    out = run_notebook("0_getting_started.ipynb", tmp_path)
    assert "Everything that DMS needs is in place" in out or "MISSING" in out
    assert "fractional luminosity f_d = 8.09e-05" in out


@requires_dms
def test_notebook_1_first_model(tmp_path):
    out = run_notebook("1_first_model.ipynb", tmp_path)
    assert "fractional luminosity f_d = 8.09e-05" in out


@requires_dms
@pytest.mark.skipif(not os.environ.get("DMS_TEST_SLOW"), reason="about two minutes; set DMS_TEST_SLOW=1")
def test_notebook_2_two_belts_observed(tmp_path):
    out = run_notebook("2_two_belts_observed.ipynb", tmp_path)
    assert "beta >= 0.44 ... 0.50" in out
    assert "f_d = 1.00e-04" in out


@requires_exe
@requires_dms
def test_notebook_3_parameter_study(tmp_path):
    """12 models with dms.batch; a second run finds all results and computes nothing."""
    for expected in ("12 models, 12 to compute", "12 models, 0 to compute"):
        out = run_notebook("3_parameter_study.ipynb", tmp_path)
        assert expected in out
    assert "ice50_G_q3.5" in out
