# Optical constants

Every file has three columns: wavelength [micron], n, k. Lines that cannot be read as three
numbers (headers, comments, a leading line count) are skipped. Fortran `D` exponents are
accepted. Use a file in the input with its path relative to this directory, e.g.
`optical_data = 'astrosilicate_draine03_862'` or
`optical_data = 'kim2019/porous_astrosilicate/as50vac50'`.

Plot files with `python tools/plot_optical_constants.py <files> -o nk.png`.

## kim2019/

Optical data of Kim et al. (2019, A&A 629, A141, "Constraining the detectability of water ice in
debris disks"). If you use them, please cite Kim et al. 2019, Reinert et al. 2015, Haßner et
al. 2016 and Potapov et al. 2019a, b (see the README.txt files in the sub-directories).

| directory | content |
|---|---|
| `astrosilicate/` | astronomical silicate (Draine 2003) |
| `ice_astrosilicate_MG/` | water ice (amorphous `a` / crystalline `c`) and astrosilicate, mixed with the Maxwell-Garnett rule; `SilXXaiceYY` = XX % silicate, YY % ice |
| `ice_mantle_astrosilicate_core/` | astrosilicate core with an ice mantle |
| `ice_platelets_astrosilicate/` | platelet-shaped ice inclusions in astrosilicate |
| `porous_astrosilicate/` | porous astrosilicate (`asXXvYY`, `asXXvacYY`: XX % silicate, YY % vacuum) |
| `porous_ices/` | porous water ice (`VacXXaiceYY`, `VacXXciceYY`) |
| `pure_ices/` | amorphous and crystalline water ice; crystalline ice for 10–260 K (temperature dependent beyond 62 micron, Reinert et al. 2015) |
| `water_ice_tables/` | source tables of the water-ice data (`.tab`, `.xml`) |
| `vacuum.txt` | vacuum (n = 1, k = 0), for porosity |

Note: `pure_ices/composite_cice_085K.txt` is identical to `composite_cice_075K.txt`, and
`ice_astrosilicate_MG/Sil30aice70` is identical to `Sil25aice75`. Please check these two files.

## Files in this directory

These are the files of the former `ri-data/` directory. The names were kept so that earlier
input files still work. Descriptions marked (?) were inferred from the file names and still need to
be confirmed.

| files | wavelength range [micron] | content |
|---|---|---|
| `astrosilicate_draine03_862`, `astrosilicate_draine03_237`, `sil100`, `draine03.txt` | 0.05–4000 (varies) | astronomical silicate (Draine 2003), different wavelength grids; `sil100` = `astrosilicate_draine03_237` |
| `silicate`, `grap-par`, `grap-per` | 0.05–2000 | silicate and graphite (parallel / perpendicular) (?) |
| `composite_aice.txt` (= `water_aice`), `composite_cice_055K.txt` (= `water_cice_55K`), `composite_cice_065K.txt`, `composite_cice_105K.txt` | 0.09–2000 | water ice (see kim2019/pure_ices) |
| `silXXciceYY`, `Sil10cice90`, `Sil50cice50_1` | 0.09–2000 / 1.1e5 | silicate–crystalline ice mixtures (?) |
| `aiXX_silYY_BR`, `ciXX_silYY_BR` | 0.09–1e4 | ice–silicate mixtures (Bruggeman rule) (?) |
| `i_Sil`, `i_ai_60K`, `i_ci_60K` | 0.09–1e4 | silicate and water ice, same wavelength grid as the `_BR` files (?) |
| `as10vac90`, `as50v50`, `as50vac50_1` | | porous astrosilicate (?) |
| `silXXforYY` | 0.1–4000 | silicate with forsterite (YY %) (?) |
| `s50mXXiYY`, `a50m50`, `dl_50_50` | 0.09–2000 | mixtures (description to be added) |
| `ch3oh.txt` | 0.05–1000 | methanol ice (?) |
