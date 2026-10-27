"""Model-atmosphere grids for DMS (&star spectrum = 'grid').

A grid directory contains ``grid_index.txt`` (lines: T_eff [K], log g [cgs], file name) and one
two-column file per model: wavelength [micron] and surface flux F_lambda [W m^-2 m^-1].
DMS interpolates bilinearly in (T_eff, log g) in log F_lambda; longward of the last tabulated
wavelength a Rayleigh-Jeans continuation is used.

Builders (internet access required for the first two):

* :func:`build_ck04_grid`   Castelli & Kurucz (2004) ATLAS9 models from the STScI archive
* :func:`build_phoenix_grid` PHOENIX-ACES-AGSS-COND-2011 (Husser et al. 2013), rebinned
* :func:`build_blackbody_grid` black bodies (tests, no download)

Command line::

    python -m dms.stellar ck04    data/stellar_models/ck04_p00 --teff 4000 10000
    python -m dms.stellar phoenix data/stellar_models/phoenix_z0 --teff 3000 7000 --logg 4.0 5.0
"""

import argparse
import io
import os
import urllib.request

import numpy as np

H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23

CK04_URLS = ("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models",
             "https://ssb.stsci.edu/cdbs/grid/ck04models")
PHOENIX_URL = "ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS"


def _write_model(directory, fname, lam_um, flam_si):
    good = np.isfinite(flam_si) & (lam_um > 0)
    np.savetxt(os.path.join(directory, fname), np.c_[lam_um[good], np.clip(flam_si[good], 0, None)],
               fmt="%.6e", header="wavelength [micron]  surface flux F_lambda [W m^-2 m^-1]")


def _write_index(directory, rows, source):
    rows = sorted(rows)
    with open(os.path.join(directory, "grid_index.txt"), "w") as f:
        f.write(f"# DMS model-atmosphere grid: {source}\n")
        f.write("# T_eff [K]   log g [cgs]   file (wavelength [micron], F_lambda [W m^-2 m^-1])\n")
        for t, g, name in rows:
            f.write(f"{t:10.1f} {g:6.2f}   {name}\n")


def _download(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def build_blackbody_grid(directory, teffs, loggs=(4.0, 4.5, 5.0), n_lambda=2000):
    """Grid of black-body spectra (for tests; log g has no effect)."""
    os.makedirs(directory, exist_ok=True)
    lam = np.logspace(np.log10(0.05), np.log10(3000.0), n_lambda)            # micron
    rows = []
    for t in teffs:
        x = H * C / (lam * 1e-6 * KB * t)
        with np.errstate(over="ignore"):
            b = 2 * H * C ** 2 / (lam * 1e-6) ** 5 / np.expm1(x)
        for g in loggs:
            name = f"bb_{int(t):05d}_{g:.2f}.dat"
            _write_model(directory, name, lam, np.pi * b)                   # F = pi B
            rows.append((float(t), float(g), name))
    _write_index(directory, rows, "black bodies")
    return directory


def build_ck04_grid(directory, metallicity="p00", teff_range=(3500, 12000), logg_range=(0.0, 5.0)):
    """Castelli & Kurucz (2004) grid, metallicity e.g. 'p00' ([M/H]=0), 'm05', 'p02'.

    Files ck{met}_{teff}.fits with columns WAVELENGTH [Angstrom] and gNN (log g = NN/10) in
    FLAM surface flux [erg s^-1 cm^-2 A^-1]; empty (all-zero) columns are skipped.
    """
    from astropy.io import fits

    os.makedirs(directory, exist_ok=True)
    teffs = [t for t in list(range(3000, 13000, 250)) + list(range(13000, 50001, 1000))
             if teff_range[0] <= t <= teff_range[1]]
    rows = []
    for t in teffs:
        data = None
        for base in CK04_URLS:
            try:
                data = _download(f"{base}/ck{metallicity}/ck{metallicity}_{t}.fits")
                break
            except Exception:                         # try the next mirror
                continue
        if data is None:
            print(f"  T = {t} K: not available, skipped")
            continue
        tab = fits.getdata(io.BytesIO(data), 1)
        lam_um = np.asarray(tab["WAVELENGTH"], float) * 1e-4
        for col in tab.columns.names:
            if not col.lower().startswith("g"):
                continue
            g = int(col[1:]) / 10.0
            if not (logg_range[0] <= g <= logg_range[1]):
                continue
            flam = np.asarray(tab[col], float)
            if not np.any(flam > 0):
                continue
            name = f"ck{metallicity}_{t:05d}_{g:.2f}.dat"
            _write_model(directory, name, lam_um, flam * 1e7)              # erg/s/cm2/A -> W/m3
            rows.append((float(t), g, name))
        print(f"  T = {t} K done")
    if not rows:
        raise RuntimeError("no models downloaded")
    _write_index(directory, rows, f"Castelli & Kurucz 2004 (ATLAS9), {metallicity}")
    return directory


def _rebin_log(lam, flux, edges):
    idx = np.searchsorted(lam, edges)
    out = np.full(len(edges) - 1, np.nan)
    cs = np.concatenate([[0.0], np.cumsum(flux)])
    for i in range(len(out)):
        a, b = idx[i], idx[i + 1]
        if b > a:
            out[i] = (cs[b] - cs[a]) / (b - a)
    return out


def build_phoenix_grid(directory, teffs, loggs, feh=0.0, n_lambda=3000):
    """PHOENIX-ACES-AGSS-COND-2011 high-resolution spectra (Husser et al. 2013), rebinned to
    n_lambda logarithmic bins between 0.05 and 5.5 micron. Fluxes: surface flux in
    erg s^-1 cm^-2 cm^-1 (converted to W m^-2 m^-1). teffs/loggs must be grid values
    (T_eff in steps of 100 K up to 7000 K, log g in steps of 0.5)."""
    from astropy.io import fits

    os.makedirs(directory, exist_ok=True)
    zdir = f"Z{'+' if feh > 0 else '-'}{abs(feh):.1f}"
    wave = fits.getdata(io.BytesIO(_download(f"{PHOENIX_URL}/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")))
    lam_um = np.asarray(wave, float) * 1e-4
    edges = np.logspace(np.log10(max(lam_um.min(), 0.05)), np.log10(lam_um.max()), n_lambda + 1)
    centre = np.sqrt(edges[1:] * edges[:-1])
    rows = []
    for t in teffs:
        for g in loggs:
            fname = f"lte{int(t):05d}-{g:.2f}{zdir[1:]}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
            try:
                flux = fits.getdata(io.BytesIO(_download(f"{PHOENIX_URL}/PHOENIX-ACES-AGSS-COND-2011/{zdir}/{fname}")))
            except Exception as exc:
                print(f"  {fname}: {exc}, skipped")
                continue
            reb = _rebin_log(lam_um, np.asarray(flux, float) * 0.1, edges)  # erg/s/cm2/cm -> W/m3
            name = f"phoenix_{int(t):05d}_{g:.2f}.dat"
            _write_model(directory, name, centre, reb)
            rows.append((float(t), float(g), name))
            print(f"  {name} done")
    if not rows:
        raise RuntimeError("no models downloaded")
    _write_index(directory, rows, f"PHOENIX-ACES-AGSS-COND-2011, [Fe/H] = {feh}")
    return directory


def main():
    p = argparse.ArgumentParser(description="Build a model-atmosphere grid for DMS")
    p.add_argument("kind", choices=["ck04", "phoenix", "blackbody"])
    p.add_argument("directory")
    p.add_argument("--teff", nargs=2, type=float, default=[3500, 10000], help="T_eff range [K]")
    p.add_argument("--logg", nargs=2, type=float, default=[3.5, 5.0], help="log g range")
    p.add_argument("--metallicity", default="p00", help="ck04: p00, m05, p02, ...")
    p.add_argument("--feh", type=float, default=0.0, help="phoenix: [Fe/H]")
    a = p.parse_args()
    if a.kind == "ck04":
        build_ck04_grid(a.directory, a.metallicity, tuple(a.teff), tuple(a.logg))
    elif a.kind == "phoenix":
        step = 100 if a.teff[1] <= 7000 else 200
        teffs = np.arange(a.teff[0], a.teff[1] + 1, step)
        loggs = np.arange(a.logg[0], a.logg[1] + 0.01, 0.5)
        build_phoenix_grid(a.directory, teffs, loggs, a.feh)
    else:
        build_blackbody_grid(a.directory, np.arange(a.teff[0], a.teff[1] + 1, 500),
                             np.arange(a.logg[0], a.logg[1] + 0.01, 0.5))


if __name__ == "__main__":
    main()
