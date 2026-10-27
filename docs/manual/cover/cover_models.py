"""Models for the cover of the DMS 2.0 manual.

Run from any directory after building the Python interface with user_density_cover.f90 in place
of src/user_density.f90 (see README.md):

    python docs/manual/cover/cover_models.py            # all models
    python docs/manual/cover/cover_models.py ring gap   # selected models

The images are written to docs/manual/cover/out/<name>.npz. NPIX (environment variable, default
401) sets the number of pixels; the cover uses 401.
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
import dms  # noqa: E402

OUT = os.path.join(HERE, "out")
NPIX = int(os.environ.get("NPIX", "401"))


def model(**groups):
    """Common settings: astronomical silicate (Draine 2003), 1-1000 micron, q = 3.5, Mie theory."""
    p = dict(run={"mode": "image", "emission": "both", "write_files": False, "screen_output": False,
                  "load_optics": True, "save_optics": True, "cache_dir": os.path.join(HERE, "cache"),
                  "data_dir": os.path.join(REPO, "data", "optical_constants")},
             star={"teff": 9500.0, "radius": 1.6, "mass": 2.2, "distance": 72.0},
             image={"n_pixels": NPIX, "n_split_max": 4},
             dust={"n_components": 1, "optical_data": "astrosilicate_draine03_862", "bulk_density": 3.5,
                   "a_min": 1.0, "a_max": 1000.0, "q": 3.5, "n_sizes": 100},
             scattering={"n_angles": 361})
    for group, values in groups.items():
        p.setdefault(group, {}).update(values)
    return dms.Model(**p)


MODELS = {
    # narrow ring in scattered and polarised light
    "ring": model(
        run={"emission": "polarisation", "optics_cache_name": "ring16"},
        disk={"density_law": "gaussian", "alpha": 0.0, "mu": 78.0, "sigma": 3.0, "r_in": 65.0, "r_out": 92.0,
              "vertical": "gaussian", "h_over_r": 0.02, "dust_mass": 1e-8,
              "inclination": 76.0, "position_angle": -62.0},
        wavelengths={"lambda": [1.6]}),
    # nearly edge-on disk
    "edge_on": model(
        run={"emission": "polarisation", "optics_cache_name": "edge16"},
        star={"teff": 8000.0, "radius": 1.5, "mass": 1.75, "distance": 19.4},
        disk={"density_law": "broken", "alpha": 5.0, "alpha_out": -3.0, "r_c": 90.0, "r_in": 20.0,
              "r_out": 250.0, "vertical": "gaussian", "h_over_r": 0.05, "dust_mass": 1e-8,
              "inclination": 87.0, "position_angle": -35.0},
        wavelengths={"lambda": [1.6]}),
    # eccentric ring: pericentre glow (24 micron) and apocentre glow (870 micron), 1.5" Gaussian PSF
    "eccentric": model(
        run={"emission": "thermal", "optics_cache_name": "ecc"},
        star={"teff": 8590.0, "radius": 1.84, "mass": 1.92, "distance": 7.7},
        disk={"density_law": "powerlaw", "alpha": 0.0, "r_in": 136.0, "r_out": 140.0, "eccentricity": 0.3,
              "omega": 0.0, "opening_angle": 2.0, "dust_mass": 5e-8,
              "inclination": 20.0, "position_angle": 25.0},
        image={"psf": "gaussian", "psf_fwhm": [1.5, 1.5]},
        wavelengths={"lambda": [24.0, 870.0]}),
    # warm inner belt and cold eccentric outer belt, three wavelengths
    "two_belts": model(
        run={"emission": "thermal", "optics_cache_name": "tb"},
        star={"teff": 8590.0, "radius": 1.84, "mass": 1.92, "distance": 7.7},
        disk={"n_belts": 2, "density_law": ["powerlaw", "gaussian"], "alpha": [0.0, 0.0], "mu": [0.0, 140.0],
              "sigma": [0.0, 10.0], "r_in": [8.0, 100.0], "r_out": [20.0, 180.0], "eccentricity": [0.0, 0.12],
              "omega": [0.0, 40.0], "belt_mass_fraction": [0.03, 0.97], "vertical": "gaussian",
              "h_over_r": 0.04, "dust_mass": 2e-7, "inclination": 55.0, "position_angle": 156.0},
        wavelengths={"lambda": [24.0, 70.0, 870.0]}),
    # broad belt with a gap (user-defined density, par(10) = 1)
    "gap": model(
        run={"emission": "thermal", "optics_cache_name": "gap"},
        star={"teff": 5850.0, "radius": 1.0, "mass": 1.0, "distance": 27.5},
        disk={"density_law": "user", "user_par": [0.0, 0.5, 80.0, 10.0, 0, 0, 0, 0, 0, 1],
              "r_in": 40.0, "r_out": 140.0, "opening_angle": 3.0, "dust_mass": 3e-8,
              "inclination": 19.0, "position_angle": 58.0},
        wavelengths={"lambda": [870.0]}),
    # ring with an azimuthal clump (the shipped example of user_density.f90, par(10) = 0)
    "clump": model(
        run={"emission": "thermal", "optics_cache_name": "clump"},
        star={"teff": 8000.0, "radius": 1.5, "mass": 1.75, "distance": 19.4},
        disk={"density_law": "user", "user_par": [1.0, 2.0, -60.0, 25.0, 0, 0, 0, 0, 0, 0],
              "r_in": 50.0, "r_out": 90.0, "opening_angle": 3.0, "dust_mass": 1e-8,
              "inclination": 30.0, "position_angle": -20.0},
        wavelengths={"lambda": [870.0]}),
}


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for name in sys.argv[1:] or list(MODELS):
        t0 = time.time()
        res = MODELS[name].run()
        images = {key: res.image(key) for key in res.images}
        np.savez_compressed(os.path.join(OUT, name + ".npz"), **images, wavelength=res.wavelength,
                            pixel=res.pixel_arcsec)
        print(f"{name}: {NPIX} pixels, {time.time() - t0:.0f} s, f_d = {res.f_d:.3e}", flush=True)
