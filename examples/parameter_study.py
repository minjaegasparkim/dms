"""DMS example 3: a parameter study with many models running in parallel (dms.batch).

How does the 3 micron water-ice band show up in the polarisation of the scattered light? The
script computes a small grid of 12 models of a nearly edge-on ring (i = 88 deg, 2.5-3.8 micron):

  3 grain compositions  astronomical silicate, and silicate with 20 % and 50 % crystalline water
                        ice (Maxwell-Garnett mixtures of Kim et al. 2019)
  2 stars               a G star and an A star; the brighter star blows out larger grains, so
                        each star has its own smallest grain size
  2 size distributions  q = 3.5 and 3.0

and plots the degree of polarisation at a scattering angle of 90 deg against wavelength, and
against the scattering angle in the ice band. (With the A star the size distribution starts at
the blow-out size, a few micron, and the sharp lower edge leaves Mie oscillations in the curves.)

Parallel computing happens on two levels:

  * inside one model, DMS uses OpenMP threads (omp_threads below = OMP_NUM_THREADS of each run);
  * dms.batch runs several models at the same time, each as its own process (workers below).

workers x omp_threads should not exceed the number of cores. (A single run from Python,
dms.Model.run, uses OpenMP alone; set OMP_NUM_THREADS before starting Python.) For many small models, many
workers with one or two threads each is fastest; for a few large models (images with many
pixels), use fewer workers with more threads. Models that differ only in the star or the
geometry share one Mie calculation (optics cache), and the runner never computes the same cache
twice at once.

The models are input files in a queue folder; every finished model keeps its own folder with the
DMS output and a result.json. If the script is stopped, running it again computes only the
models that are missing. To run a larger grid, add values to AXES; to use another computer or a
cluster, copy the queue folder there and run "python -m dms.batch run queue runs -j N" (several
runners may share one queue).

Needs the command-line program and the Python interface (make; make python). Run from the
repository root:

    python examples/parameter_study.py

About half a minute on 2 cores. Results go to parameter_study/; the figure is
parameter_study.png.
"""

import glob
import os
import sys

try:
    import dms
except ImportError:                       # not installed with pip: use the copy in this repository
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
    import dms
from dms import batch

import numpy as np
import matplotlib.pyplot as plt

OUT = "parameter_study"                   # queue/, runs/ (one folder per model) and runs/_cache/
QUEUE, RUNS = os.path.join(OUT, "queue"), os.path.join(OUT, "runs")

# ---------------------------------------------------------------------------------------------
# 1. What is common to all models
# ---------------------------------------------------------------------------------------------
base = dms.Model(
    run={"mode": "image", "emission": "polarisation", "output_format": "ascii"},
    star={"distance": 10.0},
    disk={"r_in": 80.0, "r_out": 120.0, "opening_angle": 5.0, "dust_mass": 1e-8,
          "inclination": 88.0},                  # nearly edge-on: all scattering angles are seen
    image={"n_pixels": 61, "n_split_max": 2},    # only disk-integrated quantities are used here
    wavelengths={"lambda_min": 2.5, "lambda_max": 3.8, "n_lambda": 27},
    dust={"n_components": 1, "a_min": 0.1, "a_max": 1000.0, "n_sizes": 120,   # 30 per decade
          "blowout": True},                      # grains below the blow-out size are removed
    scattering={"n_angles": 181},
)

# ---------------------------------------------------------------------------------------------
# 2. The grid: every combination of the values in AXES. A "linked" parameter changes together
#    with an axis (one value per value of the axis).
# ---------------------------------------------------------------------------------------------
COMPOSITIONS = {                                  # label: (optical data, bulk density [g/cm^3])
    "silicate": ("kim2019/astrosilicate/astrosilicate_draine2003.txt", 3.5),
    "ice20": ("kim2019/ice_astrosilicate_MG/Sil80cice20", 0.8 * 3.5 + 0.2 * 1.0),  # by volume
    "ice50": ("kim2019/ice_astrosilicate_MG/Sil50cice50", 0.5 * 3.5 + 0.5 * 1.0),
}
STARS = {"G": (5800.0, 1.0, 1.0), "A": (8590.0, 1.84, 1.92)}    # T_eff [K], R [R_sun], M [M_sun]

AXES = {
    "dust.optical_data": [c[0] for c in COMPOSITIONS.values()],
    "star.teff": [s[0] for s in STARS.values()],
    "dust.q": [3.5, 3.0],
}
LINKED = {
    "dust.optical_data": {"dust.bulk_density": [c[1] for c in COMPOSITIONS.values()]},
    "star.teff": {"star.radius": [s[1] for s in STARS.values()],
                  "star.mass": [s[2] for s in STARS.values()]},
}


def label(values):
    """Short model name, e.g. ice20_A_q3.5."""
    comp = [k for k, c in COMPOSITIONS.items() if c[0] == values["dust.optical_data"]][0]
    star = [k for k, s in STARS.items() if s[0] == values["star.teff"]][0]
    return f"{comp}_{star}_q{values['dust.q']}"


jobs = batch.grid(base, AXES, linked=LINKED, name=label)

# ---------------------------------------------------------------------------------------------
# 3. Queue the models that are not finished yet, and run them
# ---------------------------------------------------------------------------------------------
todo = [j for j in jobs if not os.path.isfile(os.path.join(RUNS, j[0], "result.json"))]
batch.write_queue(todo, QUEUE, overwrite=True)

omp_threads = 1                                   # OpenMP threads per model
workers = max(1, (os.cpu_count() or 1) // omp_threads)
print(f"{len(jobs)} models, {len(todo)} to compute: {workers} at a time, {omp_threads} thread(s) each")
batch.run(QUEUE, RUNS, workers=workers, omp_threads=omp_threads)

# ---------------------------------------------------------------------------------------------
# 4. Collect the results (one row per model) and plot
# ---------------------------------------------------------------------------------------------
rows = {r["name"]: r for r in batch.collect(RUNS)}
lam = np.array(rows[jobs[0][0]]["wavelength"])
i_band = np.abs(lam - 3.1).argmin()               # centre of the crystalline-ice band
i_cont = 0                                        # continuum at 2.5 micron

def smallest_grain(name):
    """Smallest grain radius left after blow-out, from the grain-size output file of a model."""
    f = glob.glob(os.path.join(RUNS, name, "grain_size_distribution*.dat"))[0]
    g = np.loadtxt(f, comments="#", ndmin=2)        # columns: component, a, a_low, a_high, N, ...
    return g[g[:, 4] > 0, 1].min()


print(f"\n{'model':18s} {'a_min [um]':>10s} {'f_d':>9s}   DoP(90 deg) at {lam[i_cont]:.2f} / "
      f"{lam[i_band]:.2f} micron")
for name, _ in jobs:
    r = rows[name]
    dop90 = np.array(batch.dop_at(r, 90), dtype=float)
    print(f"{name:18s} {smallest_grain(name):10.2f} {r['f_d']:9.2e}   {100 * dop90[i_cont]:5.1f} % / "
          f"{100 * dop90[i_band]:5.1f} %")

colours = {"silicate": "0.3", "ice20": "C0", "ice50": "C3"}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for k, q in enumerate([3.5, 3.0]):
    for comp in COMPOSITIONS:
        for star, ls in (("G", "-"), ("A", "--")):
            r = rows[f"{comp}_{star}_q{q}"]
            ax[k].plot(lam, 100 * np.array(batch.dop_at(r, 90), dtype=float), ls=ls,
                       color=colours[comp], label=f"{comp}, {star} star")
    ax[k].set_title(f"q = {q}: scattering angle 90 deg")
    ax[k].set_xlabel("wavelength [micron]")
    ax[k].set_ylabel("degree of polarisation [%]")
ax[0].legend(fontsize=8)

theta = np.array(rows[jobs[0][0]]["theta"])
for comp in COMPOSITIONS:
    for star, ls in (("G", "-"), ("A", "--")):
        d = np.array(rows[f"{comp}_{star}_q3.5"]["dop_theta"][i_band], dtype=float)
        ax[2].plot(theta, 100 * d, ls=ls, color=colours[comp], label=f"{comp}, {star} star")
ax[2].set_title(f"q = 3.5: {lam[i_band]:.2f} micron")
ax[2].set_xlabel("scattering angle [deg]")
ax[2].set_ylabel("degree of polarisation [%]")
ax[2].set_xlim(0, 180)

fig.tight_layout()
fig.savefig("parameter_study.png", dpi=110)
print("\nwritten: parameter_study.png")
