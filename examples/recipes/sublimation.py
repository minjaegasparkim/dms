"""Recipe: sublimation of ice (examples/sublimation_ice.nml).

Ice-silicate grains (50 % crystalline ice by volume) that become porous silicate where they are
warmer than 110 K. The sublimation radius depends on the grain size: small grains are hotter.
About 30 s.
"""
from _setup import dms, np, plt, save

res = dms.Model.from_file("examples/sublimation_ice.nml").run()
T = res.temperature                                   # T["T"][distance, component, size]
a = res.grains["a"]
r_subl = np.array([T["r_au"][T["T"][:, 0, i] >= 110.0].max() if (T["T"][:, 0, i] >= 110.0).any() else np.nan
                   for i in range(len(a))])
for i in range(0, len(a), 20):
    print(f"a = {a[i]:7.2f} micron: ice sublimates inside {r_subl[i]:5.1f} au")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for i in (0, len(a) // 3, 2 * len(a) // 3, len(a) - 1):
    ax[0].loglog(T["r_au"], T["T"][:, 0, i], label=f"{a[i]:.3g} micron")
ax[0].axhline(110, color="k", lw=0.8, ls=":")
ax[0].set_xlim(3, 150); ax[0].set_xlabel("distance [au]"); ax[0].set_ylabel("temperature of the icy grains [K]")
ax[0].legend()
ax[1].semilogx(a, r_subl, marker=".")
ax[1].set_xlabel("grain radius [micron]"); ax[1].set_ylabel("sublimation radius [au]")
save(fig, "recipe_sublimation.png")
