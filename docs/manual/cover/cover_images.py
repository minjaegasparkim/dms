"""The eight images of the cover of the DMS 2.0 manual.

Reads the model images written by cover_models.py (docs/manual/cover/out/*.npz) and writes
docs/manual/cover/images/<panel>.png: square, 1200 x 1200 pixels, black background, no text.
The cover itself (title, captions, layout) is set in LaTeX on the first page of dms_manual.tex.

Only the display is done here: every image is divided by a value close to its maximum (the
99.9th percentile of its pixels) and shown with a logarithmic or power-law stretch. The two-belt
colour composite keeps the ratios of the three wavelengths (Lupton et al. 2004, PASP 116, 133).
"""

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
IMAGES = os.path.join(HERE, "images")

BG = "#000000"
ICE = LinearSegmentedColormap.from_list("ice", ["#000000", "#0b1d3a", "#1f5d9c", "#6fb3e6", "#e8f6ff", "#ffffff"])
POL = LinearSegmentedColormap.from_list("pol", ["#000000", "#2a0b3d", "#7a1f7a", "#d24b8f", "#ffb3c7", "#ffffff"])
HEAT = plt.get_cmap("inferno")


def load(name):
    return np.load(os.path.join(OUT, name + ".npz"))


def crop_centre(img, crop):
    if crop >= 1.0:
        return img
    n = img.shape[0]
    h = int(round(n * crop / 2))
    c = n // 2
    return img[c - h:c + h + 1, c - h:c + h + 1]


def logstretch(img, decades, crop=1.0):
    img = crop_centre(np.nan_to_num(np.asarray(img, float), nan=0.0), crop)
    vmax = np.percentile(img[img > 0], 99.95)
    out = np.log10(np.clip(img, vmax * 10 ** -decades, None) / vmax) / decades + 1.0
    return np.clip(out, 0, 1)


def linstretch(img, crop=1.0, gamma=0.6):
    img = crop_centre(np.nan_to_num(np.asarray(img, float), nan=0.0), crop)
    vmax = np.percentile(img, 99.9)
    return np.clip(np.clip(img, 0, None) / vmax, 0, 1) ** gamma


def rgb_two_belts():
    """R = 870, G = 70, B = 24 micron. One logarithmic stretch of the summed intensity; the channel
    ratios (the colour) are kept (Lupton et al. 2004)."""
    d = np.nan_to_num(load("two_belts")["reem"])
    n = d.shape[1]
    y, x = np.mgrid[:n, :n] - n // 2
    outer = np.hypot(x, y) >= 40                     # pixels of the outer belt set the channel scale
    weight = (1.0, 0.5, 0.2)
    ch = [np.clip(d[i], 0, None) / np.percentile(d[i][outer], 99.9) * weight[j]
          for j, i in enumerate((2, 1, 0))]
    inten = sum(ch) / 3 + 1e-30
    decades = 2.4
    f = np.clip(np.log10(np.clip(inten / 1.5, 10 ** -decades, None)) / decades + 1, 0, 1)
    rgb = np.dstack([c * f / inten for c in ch])
    m = rgb.max(axis=2, keepdims=True)
    rgb = np.where(m > 1, rgb / np.maximum(m, 1e-30), rgb)
    white = 0.6 * f[..., None] ** 6                  # light desaturation of the brightest pixels
    return np.clip(rgb * (1 - white) + white, 0, 1)


# panel: (image, colour map); the file name is the key
PANELS = {
    "ring_sca":  (lambda: logstretch(load("ring")["sca"][0], 2.6), ICE),
    "ring_qphi": (lambda: linstretch(load("ring")["qphi"][0], gamma=0.55), POL),
    "edge_on":   (lambda: logstretch(load("edge_on")["sca"][0], 4.0, crop=0.62), ICE),
    "two_belts": (rgb_two_belts, None),
    "peri":      (lambda: linstretch(load("eccentric")["reem_obs"][0], gamma=0.8), HEAT),
    "apo":       (lambda: linstretch(load("eccentric")["reem_obs"][1], gamma=0.8), HEAT),
    "gap":       (lambda: linstretch(load("gap")["reem"][0], gamma=1.5), HEAT),
    "clump":     (lambda: linstretch(load("clump")["reem"][0], gamma=0.9), HEAT),
}


def show(ax, key):
    fn, cmap = PANELS[key]
    img = fn()
    if img.ndim == 3:
        ax.imshow(img, origin="lower", interpolation="lanczos")
    else:
        ax.imshow(img, origin="lower", cmap=cmap, vmin=0, vmax=1, interpolation="lanczos")
    ax.set_axis_off()


def write_images():
    os.makedirs(IMAGES, exist_ok=True)
    for k in PANELS:
        fig = plt.figure(figsize=(4, 4), facecolor=BG)
        ax = fig.add_axes([0, 0, 1, 1])
        show(ax, k)
        fig.savefig(os.path.join(IMAGES, f"{k}.png"), facecolor=BG, dpi=300)
        plt.close(fig)


if __name__ == "__main__":
    write_images()
