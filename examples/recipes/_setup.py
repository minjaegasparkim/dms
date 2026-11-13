"""Shared set-up of the recipes: use the DMS of this repository and run from its root.

Every recipe can be run from anywhere:  python3 examples/recipes/<recipe>.py
It writes its figure (PNG) into the current directory and prints the main numbers.
"""
import os
import sys

START = os.getcwd()                  # the figure goes here (or to DMS_RECIPE_OUT)
REPO = os.environ.get("DMS_REPO") or os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
os.chdir(REPO)                       # the example input files and cache/ are relative to the root

import dms                           # noqa: E402
import numpy as np                   # noqa: E402
import matplotlib                    # noqa: E402
if not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402

OUT = os.environ.get("DMS_RECIPE_OUT", START)


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    print("figure:", path)
