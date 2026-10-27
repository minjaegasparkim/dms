# Model-atmosphere grids

`&star spectrum = 'grid'` reads a grid directory with

- `grid_index.txt`: one line per model with `T_eff [K]   log g [cgs]   file name` (`#` starts a
  comment)
- the model files: two columns, wavelength [micron] and surface flux F_λ [W m⁻² m⁻¹]
  (F_λ = π I_λ at the stellar surface)

DMS interpolates bilinearly in (T_eff, log g) in log F_λ, on the wavelength grid of the first
model used. Longward of the last wavelength it uses a Rayleigh–Jeans continuation. The stellar
radius still comes from `&star radius`.

The grids are not stored in the repository. Build them with the Python package (internet
access needed):

```bash
# Castelli & Kurucz (2004) ATLAS9, [M/H] = 0, from the STScI archive
python -m dms.stellar ck04 data/stellar_models/ck04_p00 --teff 4000 10000 --logg 3.5 5.0

# PHOENIX-ACES-AGSS-COND-2011 (Husser et al. 2013), rebinned to 3000 points, 0.05-5.5 micron
python -m dms.stellar phoenix data/stellar_models/phoenix_z0 --teff 3000 7000 --logg 4.0 5.0

# black bodies (for tests)
python -m dms.stellar blackbody data/stellar_models/bb --teff 3000 10000
```

and use them with

```
&star  spectrum = 'grid',  grid_dir = 'data/stellar_models/ck04_p00',  teff = 5770.,  logg = 4.44,
       radius = 1.0,  mass = 1.0,  distance = 10. /
```

Please cite the model grid you use.
