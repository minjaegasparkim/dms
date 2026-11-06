#!/usr/bin/env python3
"""Check that everything DMS needs is installed, and say what to do if something is missing.

    python3 tools/check_setup.py        (or: make check)

Run it with the Python you will use for DMS (the one of your Jupyter notebook, for example):
the Python extension is compiled for one Python version and only works with that one.
Only the standard library is used, so the script runs before anything else is installed.
"""

import glob
import importlib
import importlib.util
import os
import platform
import shutil
import subprocess
import sys

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MAC = sys.platform == "darwin"
if MAC:
    INSTALL_FC = "brew install gcc   (Homebrew: https://brew.sh)"
elif sys.platform.startswith("linux"):
    INSTALL_FC = "sudo apt install gfortran   (Debian/Ubuntu; Fedora: sudo dnf install gcc-gfortran)"
else:
    INSTALL_FC = "on Windows use WSL (Ubuntu) and: sudo apt install gfortran"

rows = []          # (status, item, detail, fix)


def check(ok, item, detail="", fix="", required=True):
    status = "ok" if ok else ("MISSING" if required else "optional")
    rows.append((status, item, detail, "" if ok else fix))
    return ok


def version_of(module):
    try:
        m = importlib.import_module(module)
        return getattr(m, "__version__", "installed")
    except Exception:
        return None


def version_tuple(v):
    out = []
    for part in v.split(".")[:2]:
        digits = "".join(ch for ch in part if ch.isdigit())
        out.append(int(digits) if digits else 0)
    return tuple(out)


def main():
    pip = f"{os.path.basename(sys.executable)} -m pip install"

    # Python
    py = sys.version_info
    check(py >= (3, 9), "Python >= 3.9", f"{platform.python_version()} ({sys.executable})",
          "install Python 3.9 or newer (e.g. Anaconda or Miniforge)")

    # Fortran compiler with OpenMP
    fc = os.environ.get("FC") or "gfortran"
    path = shutil.which(fc)
    detail = ""
    if path:
        try:
            detail = subprocess.run([fc, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
        except Exception:
            detail = path
    check(bool(path), f"Fortran compiler ({fc})", detail, INSTALL_FC)
    cc = os.environ.get("CC") or "cc"
    check(bool(shutil.which(cc)), f"C compiler ({cc})", shutil.which(cc) or "",
          "macOS: xcode-select --install; Debian/Ubuntu: sudo apt install gcc "
          "(DMS 2.1 compiles one small C file)")
    check(bool(shutil.which("make")), "make", "",
          "macOS: xcode-select --install; Debian/Ubuntu: sudo apt install make")

    # NumPy and the f2py build backend
    np_v = version_of("numpy")
    check(np_v is not None and version_tuple(np_v) >= (1, 22), "NumPy >= 1.22", np_v or "",
          f"{pip} --upgrade numpy")
    if py >= (3, 12):
        meson = version_of("mesonbuild") or (shutil.which("meson") and "installed")
        ninja = shutil.which("ninja") or version_of("ninja")
        check(bool(meson), "meson (numpy.f2py with Python >= 3.12)", str(meson or ""),
              f"{pip} meson ninja")
        check(bool(ninja), "ninja (numpy.f2py with Python >= 3.12)", "installed" if ninja else "",
              f"{pip} meson ninja")

    # optional packages
    for mod, use, package in (("matplotlib", "plots and the example notebooks", "matplotlib"),
                              ("astropy", "reading FITS output", "astropy"),
                              ("pytest", "make test", "pytest"),
                              ("jupyterlab", "the example notebooks", "jupyterlab")):
        v = version_of(mod) or (version_of("notebook") if mod == "jupyterlab" else None)
        check(v is not None, mod, v or f"for {use}", f"{pip} {package}", required=False)

    # DMS itself
    exe = os.path.join(REPO, "bin", "dms")
    check(os.path.isfile(exe), "DMS program bin/dms", "", "run 'make' in the repository folder")
    ext = glob.glob(os.path.join(REPO, "python", "dms", "_dms*"))
    ext = [e for e in ext if e.endswith((".so", ".pyd"))]
    tag = f"cpython-{py.major}{py.minor}"
    ok_ext = any(tag in os.path.basename(e) for e in ext)
    detail = ", ".join(os.path.basename(e) for e in ext) or "not built"
    if ext and not ok_ext:
        fix = (f"the extension was built for another Python; run: make python PYTHON={sys.executable}")
    else:
        fix = f"run: make python PYTHON={sys.executable}"
    check(ok_ext, f"Python extension for this Python ({tag})", detail, fix)
    # another DMS installation that Python would find first (e.g. an older "pip install -e")
    spec = importlib.util.find_spec("dms")
    if spec is not None and spec.origin:
        found = os.path.normpath(os.path.dirname(spec.origin))
        here = os.path.normpath(os.path.join(REPO, "python", "dms"))
        check(found == here, "import dms finds this repository", found,
              f"another DMS is installed; remove it ({os.path.basename(sys.executable)} -m pip uninstall "
              f"dms-debris) or install this one (pip install -e {os.path.join(REPO, 'python')})")
    if ok_ext:
        sys.path.insert(0, os.path.join(REPO, "python"))
        try:
            import dms                                           # noqa: F401
            check(dms.model._dms is not None, "import dms", f"version {dms.__version__}",
                  f"rebuild: make python PYTHON={sys.executable}")
        except Exception as exc:                                 # pragma: no cover
            check(False, "import dms", str(exc), f"rebuild: make python PYTHON={sys.executable}")
    check(os.path.isfile(os.path.join(REPO, "data", "optical_constants", "astrosilicate_draine03_862")),
          "optical data (data/optical_constants)", "", "the data folder is missing: download the repository again")

    # report
    width = max(len(r[1]) for r in rows)
    print(f"DMS setup check ({platform.system()} {platform.machine()})\n")
    for status, item, detail, fix in rows:
        print(f"  {status:8s} {item:{width}s}  {detail}")
        if fix:
            print(f"  {'':8s} {'':{width}s}  -> {fix}")
    missing = [r for r in rows if r[0] == "MISSING"]
    print()
    if missing:
        print(f"{len(missing)} required item(s) missing; see the '->' lines above.")
    else:
        print("Everything that DMS needs is in place.")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
