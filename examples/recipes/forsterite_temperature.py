"""Recipe: temperature-dependent optical constants (Suto et al. 2006 forsterite, 50-295 K).

The three crystal axes of forsterite are three components, each with five temperature sets
(files data/optical_constants/suto2006/fo_suto2006_<axis>_<T>K.dat, tdep = .true.). Every grain
takes the set nearest to its own temperature ('nearest') or n, k interpolated in temperature
('linear'). Compared with the 295 K set alone, the 69 micron band of cold grains is narrower and
moves to shorter wavelength. About 1 min.
"""
from _setup import dms, np, plt, save

axes = ("A", "B", "C")
base = dms.Model(
    run={"mode": "sed", "emission": "thermal"},
    star={"teff": 8590.0, "radius": 1.84, "mass": 1.92, "distance": 20.0},
    disk={"r_in": 30.0, "r_out": 60.0, "dust_mass": 1e-8},
    wavelengths={"lambda_min": 60.0, "lambda_max": 80.0, "n_lambda": 200},
    dust={"n_components": 3, "bulk_density": [3.3] * 3, "abundance": [1.0] * 3,
          "a_min": 1.0, "a_max": 30.0, "q": 3.5, "n_sizes": 40})
cases = {
    "295 K only": {"optical_data": [f"suto2006/fo_suto2006_{x}_295K.dat" for x in axes]},
    "nearest": {"optical_data": [f"suto2006/fo_suto2006_{x}" for x in axes], "tdep": [True] * 3,
                "tdep_mode": "nearest"},
    "linear": {"optical_data": [f"suto2006/fo_suto2006_{x}" for x in axes], "tdep": [True] * 3,
               "tdep_mode": "linear", "tdep_dt": 10.0},
}
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for label, dust in cases.items():
    r = base.copy().update(dust=dust).run()
    lam, f = r.wavelength, r.sed["thermal"]
    print(f"{label:10s}  peak of the band at {lam[np.argmax(f)]:.2f} micron;  "
          f"temperature range {r.temperature['T'].min():.0f}-{r.temperature['T'].max():.0f} K; "
          f"{len(r.species['T'])} species")
    ax.plot(lam, f / np.median(f), label=label)
ax.set_xlabel("wavelength [micron]"); ax.set_ylabel("thermal flux / median")
ax.set_title("forsterite 69 micron band"); ax.legend()
save(fig, "recipe_forsterite.png")
