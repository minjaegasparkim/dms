"""The fastest recipes of examples/recipes/ (manual, Sect. 16) run and give their numbers."""
import os
import subprocess
import sys

import pytest

from conftest import REPO, requires_dms

RECIPES = os.path.join(REPO, "examples", "recipes")


def run_recipe(name, tmp_path):
    env = dict(os.environ, DMS_RECIPE_OUT=str(tmp_path), MPLBACKEND="Agg")
    out = subprocess.run([sys.executable, os.path.join(RECIPES, name)], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    return out.stdout


@requires_dms
@pytest.mark.parametrize("name, figure", [("ring_from_file.py", "recipe_ring.png"),
                                          ("eccentric_ring.py", "recipe_eccentric.png"),
                                          ("user_clump.py", "recipe_clump.png")])
def test_recipe_runs(name, figure, tmp_path):
    run_recipe(name, tmp_path)
    assert (tmp_path / figure).stat().st_size > 10000
