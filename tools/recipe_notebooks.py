"""Make the notebooks of the cookbook from the recipe scripts (same code, figures inline).

    python3 tools/recipe_notebooks.py examples/recipes
    cd examples/recipes && for n in *.ipynb; do jupyter nbconvert --to notebook --execute --inplace $n; done

Run it after changing a recipe script, so that script and notebook stay the same. The text of the
first and last cell of every notebook is in INFO below.
"""
import ast, os, re, sys
import nbformat as nbf

REC = sys.argv[1]
SETUP = '''# Setup: this notebook lives in examples/recipes/ of the repository
import os, sys
if "REPO" not in globals():
    REPO = os.environ.get("DMS_REPO") or os.path.abspath(os.path.join("..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))      # use the DMS of this repository
os.chdir(REPO)                                        # input files and cache/ are relative to the root
import dms
import numpy as np
import matplotlib.pyplot as plt
print("DMS", dms.__version__)'''

INFO = {
 "ring_from_file": ("A ring from an input file: images and SED", "16.2",
   "`examples/debris_ring.nml`", "about 10 s",
   "**Result.** f_d = 8.1e-05; 89 mJy at 70 µm and 1.8 mJy at 870 µm. The dots are the fluxes of "
   "the four images (`mode = 'both'`); the lines are the full SED of the same model, obtained by "
   "changing only the run mode and the wavelengths. From the command line, `bin/dms "
   "examples/debris_ring.nml` writes the same model to `results/` (suffix `.ring`)."),
 "eccentric_ring": ("An eccentric ring: pericentre and apocentre glow", "16.3",
   "`examples/eccentric_ring.nml`", "about 5 s",
   "**Result.** At 24 µm the pericentre side is 2.4 times brighter than the apocentre side: its "
   "dust is closer to the star and warmer (pericentre glow). At 870 µm the ratio is 0.43: the dust "
   "moves slowly near apocentre and piles up there, and the Rayleigh–Jeans emission depends only "
   "weakly on the temperature (apocentre glow)."),
 "observed_polarisation": ("Polarimetric imaging with a PSF and a coronagraph", "16.5",
   "`examples/polarisation_hg.nml`", "about 30 s",
   "**Result.** The peak of the observed Q_φ at 1.65 µm is 0.04 mJy/arcsec². For Mie polarisation "
   "instead of the Henyey–Greenstein law, remove the `&scattering` group of the input file "
   "(`phase_function = 'mie'` is the default). `examples/edge_on_disk.nml` is the same kind of "
   "observation for a disk seen at i = 88°."),
 "user_clump": ("Your own density distribution", "16.6",
   "`examples/clump_user_density.nml`, `src/user_density.f90`", "about 5 s",
   "**Result.** 9.6 mJy at 870 µm. To use another density law, write it in `user_density_fn` in "
   "`src/user_density.f90`, then run `make` and `make python` and restart the kernel."),
 "hollow_spheres": ("Non-spherical grains: hollow spheres (DHS)", "16.8",
   "model defined in Python", "about 15 s",
   "**Result.** The hollow spheres absorb more per unit mass (f_d = 5.6e-04 against 5.1e-04) and "
   "polarise very differently: the compact spheres reach −40 % (parallel to the scattering plane) "
   "at 150°, the hollow spheres +40 % near 90°."),
 "forsterite_temperature": ("Temperature-dependent optical constants", "16.9",
   "`data/optical_constants/suto2006/`", "about 35 s",
   "**Result.** The grains are at 56–89 K. With the optical constants at their temperature the "
   "69 µm band is narrower and stronger and peaks at 68.8 µm, against 69.5 µm with the 295 K set "
   "alone. `'nearest'` uses 15 data sets; `'linear'` makes 78 species with steps of at most 10 K."),
 "sublimation": ("Sublimation of ice", "16.10",
   "`examples/sublimation_ice.nml`", "about 15 s",
   "**Result.** The ice sublimates inside 52 au for 1 µm grains, which are the hottest; the "
   "sublimation radius has a minimum of 21 au near 15 µm and is 26 au for the largest grains."),
}

for stem, (title, sect, files, time, result) in INFO.items():
    src = open(os.path.join(REC, stem + ".py")).read()
    mod = ast.parse(src)
    doc = ast.get_docstring(mod)
    body_start = mod.body[1].end_lineno            # after "from _setup import ..."
    assert isinstance(mod.body[1], ast.ImportFrom)
    lines = src.splitlines()[body_start:]
    code = "\n".join(lines).strip("\n")
    code = re.sub(r'save\(fig, "[^"]+"\)', "plt.show()", code)
    blocks = [b.strip("\n") for b in re.split(r"\n\s*\n", code) if b.strip()]
    para = doc.split("\n\n", 1)
    text = para[1].replace("\n", " ") if len(para) > 1 else ""
    text = re.sub(r"\s*About [^.]*\.\s*$", "", text)
    cells = [nbf.v4.new_markdown_cell(
        f"# Recipe: {title}\n\n{text}\n\n"
        f"**Files:** {files}  \n**Run time:** {time}  \n"
        f"**Manual:** Sect. {sect} of `docs/manual/dms_manual.pdf`. The same code as a script: "
        f"`python3 examples/recipes/{stem}.py`.")]
    cells.append(nbf.v4.new_code_cell(SETUP))
    cells += [nbf.v4.new_code_cell(b) for b in blocks]
    cells.append(nbf.v4.new_markdown_cell(result))
    nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3",
         "language": "python", "name": "python3"}})
    nbf.write(nb, os.path.join(REC, stem + ".ipynb"))
    print(stem, len(cells), "cells")
