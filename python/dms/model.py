"""Python interface to the DMS Fortran code (compiled with numpy.f2py).

Example::

    import dms
    m = dms.Model.from_file("examples/debris_ring.nml")
    m.update(disk={"r_in": 60.0, "inclination": 45.0}, image={"n_pixels": 101})
    res = m.run()                          # results in memory, no files written
    res.wavelength, res.sed["thermal"]     # micron, Jy
    img = res.image("sca")                 # [n_wavelength, ny, nx] in Jy/pixel
    res.f_d, res.tau_mean

The Fortran code keeps its state in global module variables: one model runs at a time in a
Python process (OpenMP parallelises inside a run). For many models in parallel (e.g. MCMC
fitting), use separate processes (multiprocessing, MPI); each process has its own copy.
"""

import copy
import math
import os
import tempfile
import warnings

import numpy as np

from .namelist import GROUPS, read_namelist, write_namelist

try:
    from . import _dms
except ImportError as exc:                       # pragma: no cover
    _dms = None
    _import_error = exc

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_DATA = os.path.normpath(os.path.join(_PKG_DIR, "..", "..", "data", "optical_constants"))


def default_data_dir():
    """Optical-data directory: $DMS_DATA_DIR, else data/optical_constants of the repository."""
    d = os.environ.get("DMS_DATA_DIR")
    if d:
        return d
    if os.path.isdir(_REPO_DATA):
        return _REPO_DATA
    return "data/optical_constants"


class DMSError(RuntimeError):
    """Error reported by the Fortran code (invalid input, missing file, ...)."""


class SizeGridWarning(UserWarning):
    """Fewer than 25 grain-size bins per decade (see check_size_grid)."""


def check_size_grid(dust):
    """Warn if a power-law size distribution has fewer than 25 bins per decade. Each bin is
    represented by one grain size; a coarse grid leaves radial stripes in Mie Q_phi images and
    ripples with the period of the size grid in the (sub)millimetre SED. Defaults as in DMS."""
    if str(dust.get("size_distribution", "powerlaw")).lower() != "powerlaw":
        return
    n = int(dust.get("n_sizes", 100))
    a_min, a_max = float(dust.get("a_min", 1.0)), float(dust.get("a_max", 1000.0))
    if n > 1 and a_max > a_min:
        per_decade = n / math.log10(a_max / a_min)
        if per_decade < 25:
            warnings.warn(f"{per_decade:.0f} size bins per decade (n_sizes = {n}); use >= 30 per decade "
                          "for smooth Mie polarisation images and (sub)millimetre fluxes",
                          SizeGridWarning, stacklevel=3)


def _require_extension():
    if _dms is None:
        raise ImportError("the compiled DMS extension is missing; build it with 'make python' "
                          f"in the repository root ({_import_error})")


class Model:
    """A DMS model defined by the namelist groups run, star, disk, image, wavelengths, dust,
    scattering (same names and units as in the input file, see docs/input_reference.md)."""

    def __init__(self, **groups):
        unknown = set(groups) - set(GROUPS)
        if unknown:
            raise KeyError(f"unknown group(s): {', '.join(sorted(unknown))}")
        self.params = {g: dict(groups.get(g, {})) for g in GROUPS}

    @classmethod
    def from_file(cls, path):
        """Model from an existing input file (namelist)."""
        return cls(**read_namelist(path))

    def copy(self):
        return copy.deepcopy(self)

    def update(self, **groups):
        """Change parameters, e.g. ``m.update(disk={"r_in": 60}, star={"teff": 6000})``."""
        for g, entries in groups.items():
            if g not in GROUPS:
                raise KeyError(f"unknown group {g!r}")
            self.params[g].update(entries)
        return self

    def write(self, path):
        """Write the model as an input file for the command-line program."""
        write_namelist(self.params, path)

    def run(self, write_files=False, quiet=True):
        """Run DMS. Returns a Result. write_files=True also writes the usual output files."""
        _require_extension()
        params = copy.deepcopy(self.params)
        params["run"].setdefault("data_dir", default_data_dir())
        check_size_grid(params["dust"])
        fd, path = tempfile.mkstemp(suffix=".nml", prefix="dms_")
        os.close(fd)
        try:
            write_namelist(params, path)
            ierr = _dms.run_model(path, int(bool(write_files)), int(bool(quiet)))
        finally:
            os.remove(path)
        if ierr != 0:
            msg = _dms.get_error()
            msg = msg.decode() if isinstance(msg, bytes) else str(msg)
            raise DMSError(msg.strip())
        return Result._from_memory(params)

    def __repr__(self):
        used = {g: v for g, v in self.params.items() if v}
        return f"Model({used})"


class Result:
    """Results of one run (copied from the Fortran arrays).

    Attributes
    ----------
    wavelength : (n_wl,) micron
    sed : dict of (n_wl,) arrays in Jy: thermal, scattered, Q, U, polarised, star, total
    f_d, cross_section_au2, tau_mean, tau_peak, tau_radial, luminosity_lsun, teff, dust_mass_msun
    blowout_beta : (smallest, largest) blow-out threshold beta_crit; (0, 0) if blow-out is off
    pixel_arcsec, pixel_au : image sampling (images only)
    tau_profile : dict r_au, mean, pericentre, apocentre (vertical geometrical optical depth)
    grains : dict a, edges, number[n_size, n_comp], beta[n_size, n_comp]
    temperature : dict r_au, T[n_dist, n_comp, n_size]
    phase : dict angle [deg], scattered[n_wl, 180], polarised[n_wl, 180] (Jy per 1-degree bin of
            the scattering angle, whole disk)
    """

    _IMAGES = ("reem", "sca", "q", "u", "pol", "qphi", "uphi", "poldeg", "total")
    _ALIASES = {"thermal": "reem", "scattered": "sca", "polarised": "pol", "polarized": "pol",
                "degree": "poldeg"}

    @classmethod
    def _from_memory(cls, params):
        r = cls()
        r.params = params
        nwl, npix, ngr, ncomp, ndist, ntau, has_image, has_obs, complete = _dms.get_dims()
        r.n_pixels, r.has_image, r.has_obs = npix, bool(has_image), bool(has_obs)
        lam, star = _dms.get_wavelengths(nwl)
        r.wavelength = np.array(lam)
        sed = _dms.get_sed(nwl)
        r.sed = {"thermal": sed[:, 0].copy(), "scattered": sed[:, 1].copy(), "Q": sed[:, 2].copy(),
                 "U": sed[:, 3].copy(), "polarised": sed[:, 4].copy(), "star": np.array(star),
                 "total": sed[:, 5].copy()}
        v = _dms.get_scalars()
        (r.f_d, r.cross_section_au2, r.tau_mean, r.tau_peak, r.pixel_arcsec, r.pixel_au,
         r.luminosity_lsun, r.teff, r.dust_mass_msun, r.max_temperature_difference) = v[:10]
        r.tau_radial = v[11]
        r.blowout_beta = (float(v[12]), float(v[13]))      # (smallest, largest) threshold; 0 = off
        prof = _dms.get_tau_profile(ntau) if ntau > 0 else np.zeros((0, 4))
        r.tau_profile = {"r_au": prof[:, 0].copy(), "mean": prof[:, 1].copy(),
                         "pericentre": prof[:, 2].copy(), "apocentre": prof[:, 3].copy()}
        a, edges, number, beta = _dms.get_grains(ngr, ncomp)
        r.grains = {"a": np.array(a), "edges": np.array(edges), "number": np.array(number),
                    "beta": np.array(beta)}
        r_au, temp = _dms.get_temperature(ndist, ngr, ncomp)
        r.temperature = {"r_au": np.array(r_au), "T": np.transpose(temp, (2, 1, 0)).copy()}
        ang, fsca, fpol = _dms.get_phase(nwl, 180)
        r.phase = {"angle": np.array(ang), "scattered": np.array(fsca), "polarised": np.array(fpol)}
        r._images = {}
        if r.has_image:
            for name in cls._IMAGES:
                img, ok = _dms.get_image(name, npix, nwl)
                if ok:
                    r._images[name] = np.transpose(img, (2, 1, 0)).copy()      # [il, iy, ix]
                if r.has_obs:
                    img, ok = _dms.get_image(name + "_obs", npix, nwl)
                    if ok:
                        r._images[name + "_obs"] = np.transpose(img, (2, 1, 0)).copy()
            r._cubes = {}
            ang, ok = _dms.get_cube_angle(npix)
            if ok:
                r._cubes["angle"] = np.transpose(ang, (2, 1, 0)).copy()               # [iz, iy, ix]
                for name in ("reem", "sca"):
                    cube, ok = _dms.get_cube(name, npix, nwl)
                    if ok:
                        r._cubes[name] = np.transpose(cube, (3, 2, 1, 0)).copy()   # [il, iz, iy, ix]
            tau, mass = _dms.get_maps(npix)
            r.tau_map = np.array(tau).T.copy()
            r.mass_map_kg = np.array(mass).T.copy()
        return r

    def image(self, name="total", observed=False, unit="Jy/pixel"):
        """Image cube [n_wavelength, ny, nx].

        name: thermal/reem, scattered/sca, Q, U, polarised/pol, Qphi, Uphi, poldeg, total
        observed: the simulated observation (PSF / mask) if one was requested
        unit: 'Jy/pixel', 'mJy/arcsec2' or 'Jy/arcsec2' (not applied to poldeg)
        """
        key = self._ALIASES.get(name.lower(), name.lower())
        if observed:
            key += "_obs"
        if key not in self._images:
            raise KeyError(f"image {key!r} not available (available: {', '.join(self._images)})")
        img = self._images[key]
        if key.startswith("poldeg") or unit.lower() == "jy/pixel":
            return img
        factors = {"mjy/arcsec2": 1e3, "jy/arcsec2": 1.0}
        if unit.lower() not in factors:
            raise ValueError("unit must be 'Jy/pixel', 'mJy/arcsec2' or 'Jy/arcsec2'")
        scale = factors[unit.lower()] / self.pixel_arcsec ** 2
        return img * scale

    def cube(self, name="sca"):
        """3-d cube (only with image={"save_3d_cubes": True}): 'reem' or 'sca' -> flux of every
        volume element towards the observer [il, iz, iy, ix] in Jy; 'angle' -> scattering angle
        of every element [iz, iy, ix] in degrees."""
        cubes = getattr(self, "_cubes", {})
        if name not in cubes:
            raise KeyError("no 3-d cubes: run with image={'save_3d_cubes': True}")
        return cubes[name]

    def phase_curve(self, il=0, normalise=False):
        """Scattered flux of the whole disk against the scattering angle (1-degree bins).
        Returns (angle [deg], F_sca [Jy per bin], F_pol [Jy per bin]); normalise=True divides by
        the total scattered flux."""
        f, p = self.phase["scattered"][il], self.phase["polarised"][il]
        if normalise:
            tot = f.sum()
            if tot > 0:
                f, p = f / tot, p / tot
        return self.phase["angle"], f, p

    @property
    def images(self):
        """Names of the available images."""
        return list(self._images)

    def offsets(self, unit="arcsec"):
        """Pixel-centre offsets from the star along one image axis."""
        n = self.n_pixels
        k = np.arange(n) - (n - 1) / 2
        return k * (self.pixel_arcsec if unit == "arcsec" else self.pixel_au)

    def __repr__(self):
        s = f"Result(f_d={self.f_d:.3e}, tau_mean={self.tau_mean:.3e}, n_wavelength={len(self.wavelength)}"
        if self.has_image:
            s += f", images={self.n_pixels}x{self.n_pixels}"
        return s + ")"
