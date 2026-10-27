"""DMS quick start: a debris ring around a Sun-like star, computed from Python.

Run from the repository root after `make python`:

    python examples/quickstart.py

It prints a few numbers and writes quickstart.png (the SED and two images).
"""

import os
import sys

try:
    import dms
except ImportError:                       # not installed with pip: use the copy in this repository
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
    import dms

import matplotlib.pyplot as plt

# 1. The system. Names and units are those of the input files (au, degrees, M_sun, micron);
#    everything not given keeps its default (docs/input_reference.md).
model = dms.Model(
    star={"teff": 5800.0, "radius": 1.0, "mass": 1.0, "distance": 20.0},   # Sun-like star, 20 pc
    disk={"r_in": 50.0, "r_out": 100.0,             # ring from 50 to 100 au
          "dust_mass": 1e-8,                        # M_sun
          "inclination": 60.0, "position_angle": 30.0},
    dust={"n_components": 1, "optical_data": "astrosilicate_draine03_862",   # astronomical silicate
          "bulk_density": 3.5,                       # g/cm^3
          "a_min": 1.0, "a_max": 1000.0, "q": 3.5},  # grain radii [micron], n(a) ~ a^-q
)
# The grain radii are divided into n_sizes = 100 logarithmic bins (the default), each represented by
# one grain size. With fewer than ~30 bins per decade the (sub)millimetre SED ripples with the
# period of the size grid (up to ~30 per cent for 10 per decade); DMS warns below 25 per decade.

# 2. Spectral energy distribution (direct integration, no images)
sed = model.copy().update(run={"mode": "sed"},
                          wavelengths={"lambda_min": 1.0, "lambda_max": 2000.0, "n_lambda": 60}).run()
print(f"fractional luminosity f_d = {sed.f_d:.2e}")
i = sed.sed["thermal"].argmax()
print(f"thermal emission peaks at {sed.wavelength[i]:.0f} micron with {sed.sed['thermal'][i]:.3f} Jy")

# 3. Images: scattered light at 1.6 micron and thermal emission at 870 micron
img = model.copy().update(run={"mode": "image"}, wavelengths={"lambda": [1.6, 870.0]},
                          image={"n_pixels": 151}).run()
print(f"image: {img.image('total').shape[1]} x {img.image('total').shape[2]} pixels "
      f"of {img.pixel_au:.2f} au ({img.pixel_arcsec:.3f} arcsec)")

# 4. Plots
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
dms.plot.sed(sed, ax=ax[0])
top = sed.sed["total"].max()
ax[0].set_ylim(1e-5 * top, 3 * top)       # the thermal SED falls steeply below ~10 micron
dms.plot.image(img, "scattered", il=0, ax=ax[1], cmap="bone")
dms.plot.image(img, "thermal", il=1, ax=ax[2])
fig.tight_layout()
fig.savefig("quickstart.png", dpi=120)
print("written: quickstart.png")
