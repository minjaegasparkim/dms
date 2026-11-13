"""Recipe: a user-defined density - a clump in a ring (examples/clump_user_density.nml).

density_law = 'user' calls user_density_fn in src/user_density.f90 (by default a power law times
a Gaussian clump in azimuth, parameters user_par). The image is convolved with a 0.3" beam.
Edit src/user_density.f90 and run make for another density. About 15 s.
"""
from _setup import dms, np, plt, save

res = dms.Model.from_file("examples/clump_user_density.nml").run()
obs = res.image("thermal", observed=True, unit="mJy/arcsec2")
half = res.pixel_arcsec * obs.shape[-1] / 2
print(f"870 micron: {res.sed['thermal'][0] * 1e3:.2f} mJy;  peak {obs[0].max():.3f} mJy/arcsec^2")
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for j in range(2):
    p = ax[j].imshow(obs[j], origin="lower", cmap="inferno", extent=[-half, half, -half, half])
    ax[j].set_title(f"{res.wavelength[j]:g} micron, 0.3\" beam [mJy/arcsec$^2$]"); plt.colorbar(p, ax=ax[j])
    ax[j].set_xlabel("arcsec")
save(fig, "recipe_clump.png")
