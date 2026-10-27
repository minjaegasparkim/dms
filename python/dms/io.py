"""Reading the output files of the command-line program (FITS and ASCII)."""

import glob
import os

import numpy as np


def _suffix(name):
    return f".{name}" if name else ""


def read_image(output_dir, quantity="total", name="", observed=False):
    """FITS image cube [n_wavelength, ny, nx] and its header.

    quantity: reem, sca, Q, U, pol, Qphi, Uphi, poldeg, total, tau
    """
    from astropy.io import fits
    fname = os.path.join(output_dir, f"image_{quantity}{'_obs' if observed else ''}{_suffix(name)}.fits")
    with fits.open(fname) as h:
        data = np.array(h[0].data, dtype=float)
        hdr = h[0].header.copy()
    if data.ndim == 2:
        data = data[None]
    return data, hdr


def read_ascii_image(fname):
    """One ASCII image file (N rows = y, N columns = x)."""
    return np.loadtxt(fname)


def read_sed(output_dir, name=""):
    """SED file -> dict of arrays (micron, Jy)."""
    d = np.loadtxt(os.path.join(output_dir, f"SED{_suffix(name)}.dat"))
    keys = ["wavelength", "thermal", "scattered", "Q", "U", "polarised", "dust", "star", "total"]
    return {k: d[:, i] for i, k in enumerate(keys)}


def read_profile(output_dir, il=1, name=""):
    """Radial profile along the major axis for wavelength index il (1-based)."""
    return np.loadtxt(os.path.join(output_dir, f"profile_major_axis{_suffix(name)}_#{il}.dat"))


def read_optical_depth(output_dir, name=""):
    """Vertical optical depth profile: columns r [au], mean, pericentre, apocentre."""
    return np.loadtxt(os.path.join(output_dir, f"optical_depth{_suffix(name)}.dat"))


def list_outputs(output_dir, name=""):
    return sorted(os.path.basename(f) for f in glob.glob(os.path.join(output_dir, f"*{_suffix(name)}*")))
