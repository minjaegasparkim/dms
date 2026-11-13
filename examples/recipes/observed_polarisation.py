"""Recipe: a simulated polarimetric observation (examples/polarisation_hg.nml).

Near-infrared polarised light of an inclined disk with a two-component Henyey-Greenstein phase
function: diffraction-limited PSF of an 8-m telescope, a 0.1" coronagraph mask and Q_phi in
mJy/arcsec^2. About 15 s.
"""
from _setup import dms, np, plt, save

res = dms.Model.from_file("examples/polarisation_hg.nml").run()
qphi = res.image("Qphi", observed=True, unit="mJy/arcsec2")
pol = res.image("poldeg")
half = res.pixel_arcsec * qphi.shape[-1] / 2
print(f"peak Q_phi at {res.wavelength[1]:g} micron: {np.nanmax(qphi[1]):.3g} mJy/arcsec^2")
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
p = ax[0].imshow(qphi[1], origin="lower", cmap="magma", extent=[-half, half, -half, half])
ax[0].set_title(f"Q_phi observed, {res.wavelength[1]:g} micron [mJy/arcsec$^2$]"); plt.colorbar(p, ax=ax[0])
p = ax[1].imshow(100 * pol[1], origin="lower", cmap="viridis", vmin=0, vmax=60, extent=[-half, half, -half, half])
ax[1].set_title("degree of polarisation [%] (model)"); plt.colorbar(p, ax=ax[1])
for a in ax:
    a.set_xlabel("arcsec")
save(fig, "recipe_polarisation.png")
