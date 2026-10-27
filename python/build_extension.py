"""Build the Python extension dms/_dms*.so from the DMS Fortran library with numpy.f2py.

Called by ``make python`` in the repository root (which first builds build/libdms.a):

    python build_extension.py --lib ../build/libdms.a --mod ../build --omp "-fopenmp"

Requirements: NumPy (f2py), a Fortran compiler (gfortran), and for Python >= 3.12 the meson
build backend (``pip install meson ninja``).
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))

p = argparse.ArgumentParser()
p.add_argument("--lib", required=True, help="path of libdms.a")
p.add_argument("--mod", required=True, help="directory with the .mod files")
p.add_argument("--omp", default="-fopenmp", help="OpenMP flag ('' to disable)")
a = p.parse_args()

lib = os.path.abspath(a.lib)
mod = os.path.abspath(a.mod)
build_dir = os.path.join(here, "build_f2py")
shutil.rmtree(build_dir, ignore_errors=True)

env = dict(os.environ)
env["FFLAGS"] = " ".join(filter(None, [f"-I{mod}", a.omp, "-fPIC", "-O2"]))
env.setdefault("FC", "gfortran")

cmd = [sys.executable, "-m", "numpy.f2py", "-c", os.path.join(here, "dms_api.f90"),
       "-m", "_dms", "--build-dir", build_dir,
       f"-L{os.path.dirname(lib)}", "-ldms"]
if a.omp:
    cmd += ["-lgomp"] if "fopenmp" in a.omp else []
print(" ".join(cmd))
subprocess.run(cmd, check=True, cwd=here, env=env)

for f in glob.glob(os.path.join(here, "_dms*.so")) + glob.glob(os.path.join(here, "_dms*.pyd")):
    shutil.move(f, os.path.join(here, "dms", os.path.basename(f)))
print("extension placed in", os.path.join(here, "dms"))
