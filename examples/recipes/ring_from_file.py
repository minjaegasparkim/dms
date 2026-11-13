"""Recipe: a single ring from an input file - SED and images (examples/debris_ring.nml).

Reads the input file, changes nothing, runs it from Python and plots the SED and the images at
1.6 micron (scattered light) and 870 micron (thermal emission). About 20 s.
"""
from _setup import dms, np, plt, save, OUT

m = dms.Model.from_file("examples/debris_ring.nml")
res = m.run()                                    # mode = 'both': images and the SED from them
print(f"f_d = {res.f_d:.2e},  tau_mean = {res.tau_mean:.2e}")
for il, lam in enumerate(res.wavelength):
    print(f"  {lam:6.1f} micron: thermal {res.sed['thermal'][il]:.3e} Jy, scattered {res.sed['scattered'][il]:.3e} Jy")

# a full SED of the same model: change the run mode and the wavelengths only
sed = m.copy().update(run={"mode": "sed"},
                      wavelengths={"lambda_min": 0.5, "lambda_max": 2000.0, "n_lambda": 80}).run()

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
lam = sed.wavelength
for k, ls in (("thermal", "-"), ("scattered", "--"), ("star", ":")):
    ax[0].loglog(lam, sed.sed[k], ls, label=k)
ax[0].plot(res.wavelength, res.sed["thermal"] + res.sed["scattered"], "ko", ms=4, label="images (dust)")
ax[0].set_xlabel("wavelength [micron]"); ax[0].set_ylabel("flux density [Jy]"); ax[0].legend()
ax[0].set_ylim(1e-7 * sed.sed["star"].max(), 3 * sed.sed["star"].max())
half = res.pixel_arcsec * 201 / 2
for j, (name, il) in enumerate((("scattered", 0), ("thermal", 3)), start=1):
    img = res.image(name)[il]
    p = ax[j].imshow(np.log10(np.maximum(img, img.max() * 1e-4)), origin="lower", cmap="inferno",
                     extent=[-half, half, -half, half])
    ax[j].set_title(f"{name}, {lam[il]:g} micron (log Jy/pixel)")
    ax[j].set_xlabel("arcsec"); plt.colorbar(p, ax=ax[j])
save(fig, "recipe_ring.png")
