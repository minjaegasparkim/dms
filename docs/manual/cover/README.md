# Cover of the manual

The cover page is set in LaTeX at the start of `../dms_manual.tex` (title, captions and layout).
It places the eight images in `images/`, all computed with DMS 2.0. The models are described in
the section "About the cover" of the manual.

| file | content |
|---|---|
| `images/*.png` | the eight images: square, 1200 × 1200 pixels, black background, no text |
| `cover_models.py` | input of the six models (eight images) |
| `cover_images.py` | turns the model output into `images/*.png` (display stretch and colour maps) |
| `user_density_cover.f90` | user-defined density used for the gap (`user_par(10) = 1`); `user_par(10) = 0` is the shipped clump example |

To change the text or the layout of the cover, edit the `% Cover` block in `../dms_manual.tex`;
the positions are fractions of the paper width and height. Compile the manual in `docs/manual/`
with `pdflatex dms_manual.tex` (twice: the cover uses `remember picture` of TikZ).

## Reproducing the images

From the repository root:

```bash
cp docs/manual/cover/user_density_cover.f90 src/user_density.f90
make python
python docs/manual/cover/cover_models.py     # about 6 minutes on two cores
python docs/manual/cover/cover_images.py
git checkout src/user_density.f90            # restore the shipped example
make python
```

`cover_models.py` writes the model images to `out/` and keeps the optical properties in `cache/`.
Only the gap needs `user_density_cover.f90`; the other models run with the shipped code.
