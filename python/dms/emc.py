"""Effective optical constants of composite grains with the emc code (Ossenkopf 1991).

emc is an interactive Fortran program. This module writes the answers to its menus, runs it in
a temporary directory and converts the result to a DMS optical-data file (wavelength [micron],
n, k). The emc executable is taken from the argument ``emc``, else from $DMS_EMC.

Example::

    from dms import emc
    mix = [emc.Component("astrosilicate_draine03_862", 0.5, "sphere"),
           emc.Component("kim2019/pure_ices/composite_aice.txt", 0.5, "sphere")]
    emc.mix("bruggeman", mix, "data/optical_constants/mixtures/as_aice_BG_sph_f050.txt")

Theories and shapes (names as in the emc menus):

    theory               shapes allowed
    'maxwell-garnett'    'host' (exactly one), 'sphere', 'needle', 'platelet'
    'bruggeman'          'sphere', 'needle', 'platelet'
    'uni-mg'             'host' (exactly one), 'sphere', 'needle', 'platelet',
                         'coated' (core; give mantle=...), 'coating' (the mantle material),
                         'shape' (ellipsoid with depolarisation factor L)
    'uni-bg'             as 'uni-mg' without 'host'

Vacuum is emc's built-in material 7: give the file name 'vacuum'. Input files are DMS optical
data files (wavelength in micron); emc works in cm, the conversion is done here. The output
uses the wavelength grid of the first non-vacuum component, restricted to the range covered by
all components.
"""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Optional

import numpy as np

THEORIES = {"maxwell-garnett": 1, "bruggeman": 2, "uni-mg": 3, "uni-bg": 4}
_ALIASES = {"mg": "maxwell-garnett", "bg": "bruggeman", "unimg": "uni-mg", "unibg": "uni-bg",
            "universal-mg": "uni-mg", "universal-bg": "uni-bg"}

# emc menu numbers of the shapes, per theory
_SHAPES = {
    1: {"host": 1, "sphere": 2, "needle": 3, "platelet": 4},
    2: {"sphere": 1, "needle": 2, "platelet": 3},
    3: {"sphere": 1, "needle": 2, "platelet": 3, "coated": 4, "coating": 5, "shape": 6, "host": 7},
    4: {"sphere": 1, "needle": 2, "platelet": 3, "coated": 4, "coating": 5, "shape": 6},
}


class EMCError(RuntimeError):
    """emc did not produce a result (wrong input, wavelength range, volume fractions ...)."""


@dataclass
class Component:
    """One material of the mixture.

    file:     DMS optical-data file (relative to data_dir, or absolute), or 'vacuum'
    fraction: volume fraction (all fractions must add up to 1)
    shape:    see the module docstring
    mantle:   for shape='coated': index (0-based) of the component that forms the mantle
    L:        for shape='shape': depolarisation factor (0..1)
    """
    file: str
    fraction: float
    shape: str = "sphere"
    mantle: Optional[int] = None
    L: Optional[float] = None


def _theory_number(theory):
    t = theory.lower()
    t = _ALIASES.get(t, t)
    if t not in THEORIES:
        raise ValueError(f"unknown theory {theory!r}; use one of {', '.join(THEORIES)}")
    return THEORIES[t]


def _read_nk(path):
    """Wavelength [micron], n, k from a DMS optical-data file (comment lines allowed)."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#!%":
                continue
            try:
                vals = [float(v) for v in s.replace(",", " ").split()[:3]]
            except ValueError:
                continue
            if len(vals) == 3:
                rows.append(vals)
    a = np.array(rows)
    if a.ndim != 2 or len(a) < 2:
        raise EMCError(f"no data in {path}")
    return a[np.argsort(a[:, 0])]


def mix(theory, components, output, data_dir=None, emc=None, header=None, keep_log=False, grid_from=None):
    """Compute effective n, k with emc and write them to ``output`` (DMS format).

    grid_from: index of the component whose wavelength grid is used for the result (default: the
    first component that is not vacuum); the result covers the range shared by all components.
    Returns the array [wavelength (micron), n, k].
    """
    tn = _theory_number(theory)
    exe = emc or os.environ.get("DMS_EMC")
    if not exe or not os.path.isfile(exe):
        raise EMCError("emc executable not found: give emc=... or set $DMS_EMC")
    exe = os.path.abspath(exe)
    if data_dir is None:
        from .model import default_data_dir
        data_dir = default_data_dir()

    comps = [c if isinstance(c, Component) else Component(**c) for c in components]
    if not 1 <= len(comps) <= 6:
        raise ValueError("emc mixes 1 to 6 components")
    ftot = sum(c.fraction for c in comps)
    if abs(ftot - 1.0) > 1e-6:
        raise ValueError(f"volume fractions add up to {ftot}, not 1")
    shapes = _SHAPES[tn]
    for i, c in enumerate(comps):
        if c.shape not in shapes:
            raise ValueError(f"shape {c.shape!r} not available for theory {theory!r} "
                             f"(allowed: {', '.join(shapes)})")
        if c.shape == "coated":
            if c.mantle is None or not (0 <= c.mantle < len(comps)) or c.mantle == i:
                raise ValueError(f"component {i}: 'coated' needs mantle=<index of another component>")
            if comps[c.mantle].shape != "coating":
                raise ValueError(f"component {c.mantle} is the mantle of component {i}: give it shape='coating'")
        if c.shape == "shape" and (c.L is None or not 0.0 <= c.L <= 1.0):
            raise ValueError(f"component {i}: shape='shape' needs 0 <= L <= 1")
    if "host" in shapes and tn in (1, 3):
        nhost = sum(c.shape == "host" for c in comps)
        if nhost != 1:
            raise ValueError("Maxwell-Garnett theories need exactly one component with shape='host'")
    if all(c.file.lower() == "vacuum" for c in comps):
        raise ValueError("at least one component must not be vacuum")
    if sum(c.file.lower() == "vacuum" for c in comps) > 1:
        raise ValueError("emc has a single vacuum material: use 'vacuum' for one component only "
                         "(a host of fraction 0 can be any other material file)")

    with tempfile.TemporaryDirectory(prefix="emc_") as tmp:
        seq, matno, nmat, src = [], [], 0, []
        for c in comps:
            if c.file.lower() == "vacuum":
                matno.append(7)
                continue
            nmat += 1
            path = c.file if os.path.isabs(c.file) else os.path.join(data_dir, c.file)
            d = _read_nk(path)
            src.append((os.path.abspath(path), d[0, 0], d[-1, 0]))
            np.savetxt(os.path.join(tmp, f"m{nmat}.lnk"),
                       np.column_stack([d[:, 0] * 1e-4, d[:, 1], d[:, 2]]), fmt="%.8e")
            seq += ["1", str(nmat), f"m{nmat}.lnk"]
            matno.append(nmat)
        if nmat >= 6:
            raise ValueError("at most 5 non-vacuum components (emc keeps the result in slot 6)")
        res = nmat + 1
        seq += ["3", str(res)]
        for m in matno:
            seq += ["1", str(m)]
        for m, c in zip(matno, comps):
            seq += ["3", str(m), repr(float(c.fraction))]
        if grid_from is None:
            first = next(m for m in matno if m != 7)
        else:
            first = matno[grid_from]
            if first == 7:
                raise ValueError("grid_from must not point to the vacuum component")
        seq += ["4", "1", str(first)]
        seq += ["5", str(tn)]
        for i, c in enumerate(comps):
            seq.append(str(shapes[c.shape]))
            if c.shape == "coated":
                seq.append(str(matno[c.mantle]))
            if c.shape == "shape":
                seq.append(repr(float(c.L)))
        seq += ["n"]                      # no correction of the shapes
        seq += ["2", str(res), "out.lnk", "5", "y"]
        stdin = "\n".join(seq) + "\n"
        try:
            p = subprocess.run([exe], input=stdin, cwd=tmp, capture_output=True, text=True,
                               timeout=600)
        except subprocess.TimeoutExpired:
            raise EMCError("emc did not finish within 10 minutes")
        log = p.stdout + p.stderr
        if keep_log:
            with open(str(output) + ".emc.log", "w") as f:
                f.write("--- input ---\n" + stdin + "--- output ---\n" + log)
        for msg in ("error occured", "There is no access", "not present", "cannot", "Type 0-", "Type 1-"):
            if msg in log:
                raise EMCError(f"emc reported a problem ({msg!r}); run with keep_log=True to see its output")
        outp = os.path.join(tmp, "out.lnk")
        if p.returncode != 0 or not os.path.isfile(outp):
            raise EMCError(f"emc failed (exit code {p.returncode}); run with keep_log=True")
        out = np.loadtxt(outp)
    out = out[np.argsort(out[:, 0])]
    lam = out[:, 0] * 1e4
    res = np.column_stack([lam, out[:, 1], out[:, 2]])
    if not np.all(np.isfinite(res)) or np.any(res[:, 1] <= 0) or np.any(res[:, 2] < 0):
        raise EMCError("emc returned non-finite or unphysical values (n <= 0 or k < 0)")

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    lines = ["# effective optical constants computed with emc (Ossenkopf 1991)",
             f"# theory: {theory}"]
    for i, c in enumerate(comps):
        extra = f", mantle = component {c.mantle}" if c.shape == "coated" else ""
        extra += f", L = {c.L}" if c.shape == "shape" else ""
        lines.append(f"# component {i}: {c.file}, volume fraction {c.fraction}, {c.shape}{extra}")
    if header:
        lines += ["# " + h for h in header.splitlines()]
    lines.append("# wavelength [micron]   n   k")
    tmpout = str(output) + ".tmp%d" % os.getpid()
    with open(tmpout, "w") as f:
        f.write("\n".join(lines) + "\n")
        np.savetxt(f, res, fmt="%.6e")
    os.replace(tmpout, output)            # atomic: parallel jobs never see a half-written file
    return res


def bulk_density(components, densities):
    """Volume-weighted bulk density of the mixture (vacuum has density 0).

    densities: dict file -> g/cm^3 (vacuum may be omitted)."""
    rho = 0.0
    for c in components:
        c = c if isinstance(c, Component) else Component(**c)
        if c.file.lower() != "vacuum":
            rho += c.fraction * densities[c.file]
    return rho
