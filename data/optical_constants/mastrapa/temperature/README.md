# Water ice at several temperatures

Files provided by the author (2026-10-23), unchanged. Columns: wavelength [micron], n, k;
0.091–10⁴ micron, 2298 points on the same wavelength grid in every file (0.005 micron steps in
2.7–3.4 micron).

| files | phase | temperatures [K] |
|---|---|---|
| `i_ai_<T>K` | amorphous | 15, 25, 40, 60, 80, 100, 120 |
| `i_ci_<T>K` | crystalline | 20, 30, 40, 60, 80, 100, 120, 140 |

The temperatures are those of the Mastrapa et al. (2008, 2009) measurements (the 50 K data are in
`../aice_50K_hybrid.txt`, `../cice_50K_hybrid.txt`). Outside the measured range the data are joined
with those of Kim et al. (2019), whose temperature dependence lies at far-infrared and
(sub)millimetre wavelengths (Mishima et al. 1983; Reinert et al. 2015; Häßner et al. 2018). `i_ai_60K` and `i_ci_60K` are identical to the files of the same
name in `data/optical_constants/`.

Checks (2026-10-23)

- The data depend on temperature only at 0.55–16 micron (amorphous) and 1.2–16 micron
  (crystalline); in 1.1–2.5 micron the change is below 0.01 in n and 0.0003 in k.
- 3 micron band: amorphous ice, the k peak grows with T (0.46 at 15 K to 0.72 at 120 K) and moves
  from 3.035 to 3.080 micron; crystalline ice, the k peak falls with T (1.22 at 20 K to 0.87 at
  140 K) at 3.10–3.115 micron.
- Crystalline files: n steps from 1.320 to 1.309 between 1.285 and 1.295 micron (a join of two
  data sets; the 50 K hybrid file has n = 1.309 there without a step). Amorphous files: smooth.
- Against the 50 K files (1.3–5 micron, data interpolated to 0.005 micron steps): crystalline
  50 K lies within the 40–60 K range to 0.002 in n and 0.003 in k. Amorphous 50 K lies outside
  the 40–60 K range in parts of 2.57–4.9 micron, by up to 0.032 in n (3.1 micron) and 0.027 in k
  (3.24 micron). Below 1.29 micron the crystalline files differ from the 50 K file by the step
  above (0.014 in n).
