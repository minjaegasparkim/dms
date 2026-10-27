#!/usr/bin/env python3
"""Quick-look figure of a finished DMS run (FITS output): images at all wavelengths, the SED and
the vertical optical depth.

    python tools/quicklook.py results ring -o ring.png        (output directory, output_name)

(Replaces the scripts of the former visualization/ folder, which read the old ASCII outputs.)
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from dms import io  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("output_dir")
    p.add_argument("name", nargs="?", default="")
    p.add_argument("-q", "--quantity", default="total", help="reem, sca, pol, Qphi, total, ... (default total)")
    p.add_argument("--observed", action="store_true", help="use the *_obs images")
    p.add_argument("-o", "--output", default="dms_quicklook.png")
    a = p.parse_args()

    img, hdr = io.read_image(a.output_dir, a.quantity, a.name, a.observed)
    nwl = img.shape[0]
    lam = [hdr.get(f"WAVE{i + 1:04d}", np.nan) for i in range(nwl)]
    pix = hdr["CDELT1"]
    n = img.shape[-1]
    ext = (np.array([-0.5, n - 0.5]) - (n - 1) / 2) * pix
    extent = [ext[0], ext[1], ext[0], ext[1]]

    has_sed = os.path.exists(os.path.join(a.output_dir, f"SED{'.' + a.name if a.name else ''}.dat"))
    ncol = nwl + 1 + int(has_sed)
    fig, ax = plt.subplots(1, ncol, figsize=(3.6 * ncol, 3.4))
    for i in range(nwl):
        d = img[i]
        vmax = np.nanmax(d)
        ax[i].imshow(np.log10(np.clip(d, vmax * 1e-4, None)) if vmax > 0 else d, origin="lower",
                     extent=extent, cmap="inferno")
        ax[i].set_title(f"{a.quantity} {lam[i]:g} um", fontsize=9)
        ax[i].set_xlabel("arcsec")
    k = nwl
    if has_sed:
        s = io.read_sed(a.output_dir, a.name)
        for key in ("thermal", "scattered", "star", "total"):
            y = s[key]
            ax[k].loglog(s["wavelength"], np.where(y > 0, y, np.nan), ".-", label=key)
        top = np.nanmax(s["total"])
        ax[k].set_ylim(top * 1e-8, top * 3)
        ax[k].set_xlabel("wavelength [micron]")
        ax[k].set_ylabel("Jy")
        ax[k].legend(fontsize=7)
        k += 1
    t = io.read_optical_depth(a.output_dir, a.name)
    ax[k].plot(t[:, 0], t[:, 1], label="azimuthal mean")
    if np.any(t[:, 2] != t[:, 1]):
        ax[k].plot(t[:, 0], t[:, 2], label="pericentre")
        ax[k].plot(t[:, 0], t[:, 3], label="apocentre")
    ax[k].set_xlabel("r [au]")
    ax[k].set_ylabel("vertical tau")
    ax[k].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(a.output, dpi=110)
    print("written", a.output)


if __name__ == "__main__":
    main()
