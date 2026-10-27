"""Quick-look plots (matplotlib) for Result objects or output files."""

import numpy as np


def image(res, name="total", il=0, observed=False, unit="Jy/pixel", log=True, ax=None, **kw):
    """Show one image of a Result (offsets in arcsec, star at (0, 0))."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    img = res.image(name, observed=observed, unit=unit)[il]
    ext = res.offsets("arcsec")
    h = 0.5 * res.pixel_arcsec
    extent = [ext[0] - h, ext[-1] + h, ext[0] - h, ext[-1] + h]
    if log:
        vmax = np.nanmax(img)
        img = np.log10(np.clip(img, vmax * 1e-4, None)) if vmax > 0 else img
    im = ax.imshow(img, origin="lower", extent=extent, cmap=kw.pop("cmap", "inferno"), **kw)
    ax.set_xlabel("offset x [arcsec]")
    ax.set_ylabel("offset y [arcsec]")
    ax.set_title(f"{name} {res.wavelength[il]:g} micron" + (" (observed)" if observed else ""))
    plt.colorbar(im, ax=ax, label=("log " if log else "") + (unit if name != "poldeg" else "P"))
    return ax


def sed(res, ax=None, components=("thermal", "scattered", "star", "total")):
    """lambda F_lambda-like plot of the SED in Jy."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    for k in components:
        y = res.sed[k]
        ax.loglog(res.wavelength, np.where(y > 0, y, np.nan), label=k, marker=".")
    ax.set_xlabel("wavelength [micron]")
    ax.set_ylabel("flux density [Jy]")
    ax.legend()
    return ax


def temperature(res, component=0, ax=None, sizes=None):
    """Grain temperature as a function of distance for some grain sizes."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    a = res.grains["a"]
    idx = sizes if sizes is not None else np.unique(np.linspace(0, len(a) - 1, 5).astype(int))
    for i in idx:
        ax.plot(res.temperature["r_au"], res.temperature["T"][:, component, i], label=f"a = {a[i]:.3g} micron")
    ax.set_xlabel("distance [au]")
    ax.set_ylabel("T [K]")
    ax.legend()
    return ax
