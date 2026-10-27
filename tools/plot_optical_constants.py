#!/usr/bin/env python3
"""Plot n and k of one or more optical-data files (columns: wavelength [micron], n, k).

    python tools/plot_optical_constants.py data/optical_constants/astrosilicate_draine03_862 \
        data/optical_constants/kim2019/pure_ices/composite_cice_065K.txt -o nk.png

(Replaces the n_k.py scripts that were stored next to the data files.)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(fname):
    rows = []
    with open(fname) as f:
        for line in f:
            try:
                v = [float(x.replace("D", "E").replace("d", "e")) for x in line.split()[:3]]
            except ValueError:
                continue
            if len(v) == 3:
                rows.append(v)
    return np.array(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("-o", "--output", default="optical_constants.png")
    a = p.parse_args()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for f in a.files:
        d = load(f)
        lab = os.path.basename(f)
        ax[0].semilogx(d[:, 0], d[:, 1], label=lab)
        ax[1].loglog(d[:, 0], d[:, 2], label=lab)
    for x, yl in zip(ax, ["n", "k"]):
        x.set_xlabel("wavelength [micron]")
        x.set_ylabel(yl)
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.output, dpi=120)
    print("written", a.output)


if __name__ == "__main__":
    main()
