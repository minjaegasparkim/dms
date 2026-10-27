# Refractory materials from optool

Optical constants (wavelength [micron], n, k) taken from the `lnk_data` folder of optool
(Dominik, Min & Tazaki 2021, Astrophysics Source Code Library, ascl:2104.010; MIT licence, see
`LICENSE_optool.txt`; repository github.com/cdominik/optool, commit a38fd6a, 2026-10-12).
Only the format was changed: the first line of the optool files (number of points and bulk
density) was moved into a comment.

| file | material | bulk density [g/cm³] | wavelengths [µm] | reference |
|---|---|---|---|---|
| `pyr-mg70-Dorschner1995.txt` | amorphous pyroxene Mg0.7Fe0.3SiO3 | 3.01 | 0.2–500 | Dorschner et al. (1995), A&A 300, 503 |
| `c-z-Zubko1996.txt` | amorphous carbon | 1.80 | 0.05–10⁴ | Zubko et al. (1996), MNRAS 282, 1321 |
| `c-org-Henning1996.txt` | organics (CHON) | 1.5 | 0.1–10⁴ (cut here; the optool file reaches 10⁵) | Henning & Stognienko (1996), A&A 311, 291; originally Pollack et al. (1994) |
| `diana_pyr70_c_BG.txt` | pyroxene + amorphous carbon, mass ratio 0.87 : 0.13 (Bruggeman) | 2.768 | 0.2–500 | as the DIANA standard dust (Woitke et al. 2016, A&A 586, A103), without porosity |

Notes

- Zubko et al. (1996) measured several carbon samples (ACAR, ACH2, BE). The optool file does
  not say which one it contains, and its table lists the state as "amorph?". Check before
  citing a sample.
- Organics: the optool file gives 1.5 g/cm³, the optool user guide table 1.4 g/cm³.
- Pyroxene has 20 points and organics 24 points between 1 and 5 µm. Both curves are smooth
  there, and emc interpolates them to the grid of the ice data.
- `diana_pyr70_c_BG.txt` was computed with emc (Bruggeman, spheres): volume fractions
  0.8001 (pyroxene) and 0.1999 (carbon), using the densities above.
