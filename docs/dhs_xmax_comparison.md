# DHS for large grains: `dhs_xmax` (DMS 2.1.1)

Measurements behind the default `dhs_xmax = 1000`. Machine: 2 threads (gfortran 13, `-O3
-march=native`). The test model is a debris ring around an A star (27.75–32.25 au, 19.44 pc),
astronomical silicate + forsterite axes A, B, C (Suto et al. 2006, 50–295 K, `'nearest'`,
15 sets), `grain_shape = 'dhs'`, `dhs_fmax = 0.8`, 10 points in f, 0.264–a_max micron, eight
observing wavelengths (10–100 micron), thermal emission only.

## Where the time goes

| step | time (a_max = 100 micron, 30 sizes, full DHS) |
|---|---|
| optical-data grid (2792 wavelengths x 30 sizes x 16 species x 10 f) | 218 s |
| observing wavelengths | 0.1 s |
| temperatures | 0.9 s |
| rest | < 1 s |

One coated-sphere call (hollow sphere, f = 0.45) against one Mie call (miex) of the same x:

| x | 1 | 10 | 100 | 1e3 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|
| coated sphere [ms] | 0.0018 | 0.0038 | 0.021 | 0.20 | 2.0 | 20 |
| Mie [ms] | 0.0004 | 0.0011 | 0.0056 | 0.048 | 0.47 | 4.4 |

Both grow in proportion to x. With grains up to 1 mm and optical data from 0.1 micron, the
grains with x > 1000 are 11 % of the (wavelength, size) points and 90 % of the time. The loop over
wavelengths, sizes and species runs in parallel (OpenMP); the loop over f inside one grain is
serial, which is not limiting (millions of parallel iterations).

## Approximation

A grain with x > `dhs_xmax` is computed at x = `dhs_xmax` with k' = k x / `dhs_xmax` (same k x),
limited to k' `dhs_xmax` = 30 (opaque); if k `dhs_xmax` >= 30, k is kept. Replacing large DHS grains
by volume-equivalent compact spheres was rejected: Qabs of the sphere is 36 % lower for absorbing
material (the hollow spheres have a larger cross section per volume, on average (1-f)^(-2/3) =
1.56 for f_max = 0.8) and up to 60 % higher for transparent forsterite.

Efficiencies against the full DHS, x = 3e3–7e4 (largest relative difference):

| material | x_max = 300 | 1000 | 3000 |
|---|---|---|---|
| astronomical silicate: Qabs, Qpr, Qsca | 1.9 %, 1.7 %, 1.9 % | 0.65 %, 0.55 %, 0.89 % | 0.22 %, 0.19 %, 0.42 % |
| forsterite C, 295 K | 7.3 %, 2.4 %, 3.5 % | 4.3 %, 1.3 %, 1.9 % | 1.5 %, 0.45 %, 0.60 % |
| forsterite A, 50 K | 12 %, 2.3 %, 4.4 % | 14 %, 1.2 %, 1.9 % | 1.5 %, 0.48 %, 0.82 % |

The large differences for forsterite are at 0.3–1 micron, where k = 5e-8–1e-6 and Qabs < 1e-3. The
full DHS itself changes by up to ~10 % there between 10 and 40 points in f.

## Model results against the full DHS

Relative difference of the thermal flux, and of the grain temperatures (all components, sizes and
distances):

| a_max, sizes | dhs_xmax | 10 | 33.5 | 44.0 | 49.3 | 62.0 | 69.0 | 69.2 | 100 micron | max dT/T | time |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 100, 30 | full | | | | | | | | | | 191 s |
| 100, 30 | 3000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 172 s |
| 100, 30 | 1000 | 0.0004 | 0.0001 | 0.0001 | 0.0001 | 0.0001 | 0.0001 | 0.0001 | 0.0001 | 0.0007 | 112 s |
| 100, 30 | 300 | 0.0088 | 0.0020 | 0.0016 | 0.0015 | 0.0014 | 0.0013 | 0.0013 | 0.0014 | 0.0029 | 72 s |
| 1000, 30 | full | | | | | | | | | | 1204 s |
| 1000, 30 | 1000 | 0.0036 | 0.0005 | 0.0004 | 0.0004 | 0.0004 | 0.0003 | 0.0003 | 0.0003 | 0.0014 | 247 s |
| 1000, 30 | 300 | 0.0192 | 0.0036 | 0.0028 | 0.0026 | 0.0024 | 0.0023 | 0.0023 | 0.0023 | 0.0043 | 119 s |

All differences are positive (slightly warmer grains). The full DHS of version 2.1.0 took 220 s for
the first model (the coated-sphere routine is now ~15 % faster). `grain_shape = 'sphere'` gives the
same results bit for bit (one thread).

## Test model with 100 sizes, a_max = 1000 micron

| | time (2 threads) |
|---|---|
| version 2.1.0 (full DHS) | ~71 min (estimated: 760 s x 1204 s / 247 s x 1.15) |
| 2.1.1, `dhs_xmax = 1000` | 760 s |
| 2.1.1 with 3233 observing wavelengths (24–235 micron), first run | 982 s (optical-data grid 750 s, observing wavelengths 167 s, rest 66 s) |
| same, optics loaded from the cache | 67 s (emission kernels 61 s) |
