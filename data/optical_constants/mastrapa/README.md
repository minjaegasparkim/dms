# Water ice at 50 K (Mastrapa et al.)

| file | content |
|---|---|
| `aice_50K_hybrid.txt` | amorphous water ice, 50 K |
| `cice_50K_hybrid.txt` | crystalline water ice, 50 K |

Columns: wavelength [micron], n, k.

- 1.111–22.23 micron (amorphous) and 1.111–15.00 micron (crystalline): laboratory data of
  Mastrapa et al. at 50 K, measured with a constant wavenumber step of about 0.965 cm^-1. Points
  are kept with their original values; only some points are left out, so that the spacing is
  Δλ/λ ≥ 2·10⁻⁴ in 2.7–3.5 micron, ≥ 6·10⁻⁴ in 1.9–5.5 micron and ≥ 5·10⁻³ elsewhere
  (linear interpolation of the thinned data differs from the full data by less than 0.0003 in
  n and k at 2–5 micron).
- Outside these ranges: `kim2019/pure_ices/composite_aice.txt` and
  `kim2019/pure_ices/composite_cice_055K.txt` (Kim et al. 2019), which DMS needs for the
  temperatures (stellar and thermal wavelengths). The two data sets do not join smoothly: at
  22.2 micron (amorphous) n = 1.75 / 1.47 and k = 0.137 / 0.016; at 15.0 micron (crystalline)
  n = 1.64 / 1.52 and k = 0.120 / 0.305 (Mastrapa / Kim et al. 2019).

The first column is the wavelength in micron (as all DMS optical data); emc reads wavelengths in
cm, `dms.emc` converts.
