"""Figure 5: higher-order dispersion and entanglement structure.
(a) repetition-rate shift vs D3 scaling with crystal-stability boundary,
(b) drop-port squeezing spectra for D3 scalings,
(c) supermode weight asymmetry vs D3 scaling,
(d) drop-port logarithmic-negativity matrix of the odd-mode lattice.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import SERIES, OI, FULLW, panel_label

d3q = np.load("q_d3.npz")
d3f = np.load("lle_d3family.npz")
ent = np.load("q_entanglement.npz")
dat = np.load("lle_family.npz")
kappa = float(dat["kappa"])
D3_design = float(dat["D"][2])

fig, axs = plt.subplots(2, 2, figsize=(FULLW, 3.65))

scales = np.array([float(k) for k in d3f["keys"]])
vs = np.array(d3f["v"])

# (a) rep-rate shift vs D3
ax = axs[0, 0]
d3_khz = scales * D3_design / 2 / np.pi / 1e3
ax.plot(d3_khz, vs * (kappa / 2) / 2 / np.pi / 1e6, "o-", ms=3.5,
        color=SERIES[2])
bnd = float(d3f["boundary_scale"]) * D3_design / 2 / np.pi / 1e3
ax.axvspan(bnd, 1.15 * bnd, color=OI["vermilion"], alpha=0.18, lw=0)
ax.axvline(bnd, color=OI["vermilion"], lw=0.8, ls="--")
ax.text(bnd * 1.045, -2.2, "crystal lost", rotation=90, fontsize=6.4,
        ha="center", va="center", color=OI["vermilion"])
ax.axvline(D3_design / 2 / np.pi / 1e3, color="0.5", lw=0.7, ls=":")
ax.annotate("FEM value", xy=(D3_design / 2 / np.pi / 1e3, -3.55),
            xytext=(-120, -4.15), fontsize=6.4,
            arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
ax.set_xlabel(r"$|D_3|/2\pi$ (kHz)")
ax.set_ylabel(r"$\delta f_{\rm rep}$ (MHz)")
ax.grid(True)
panel_label(ax, "(a)")

# (b) spectra vs D3: curves are indistinguishable at full scale, so show
# one representative curve and resolve the family in a zoom inset
ax = axs[0, 1]
oms = d3q["oms"]
ax.plot(oms, 10 * np.log10(d3q["smin_0.0"]), lw=1.2, color="0.25")
ax.axhline(0, color="0.5", lw=0.6)
ax.set_xlim(0, 40)
ax.set_ylim(-8.6, 0.6)
ax.set_xlabel(r"sideband frequency $\omega/\kappa$")
ax.set_ylabel("squeezing (dB)")
ax.annotate("all $D_3$ scales\n(0--3$\\times$, overlapping)",
            xy=(8.2, -5.1), xytext=(2.3, -7.9), fontsize=6.0,
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.4"))
# inset: zoom on the spectral minimum, 5 curves resolved
axi = ax.inset_axes([0.52, 0.13, 0.45, 0.52])
for j, s3 in enumerate(scales):
    axi.plot(oms, 10 * np.log10(d3q[f"smin_{s3}"]), lw=0.9,
             color=SERIES[j % 6], label=f"{s3:g}$\\times$")
axi.set_xlim(0, 4)
axi.set_ylim(-7.302, -7.243)
axi.set_yticks([-7.30, -7.25])
axi.set_yticklabels([])
axi.tick_params(labelsize=5, pad=1.5, length=2)
axi.set_title("zoom: spread $<$ 0.01 dB", fontsize=5.4, pad=1.5)
axi.legend(fontsize=4.6, ncol=2, loc="upper right", borderpad=0.25,
           columnspacing=0.6, handletextpad=0.35, handlelength=1.1,
           title="$D_3$ scale", title_fontsize=4.8)
ax.indicate_inset_zoom(axi, edgecolor="0.5", lw=0.6)
panel_label(ax, "(b)")

# (c) supermode asymmetry
ax = axs[1, 0]
odd_mu = np.array([m for m in range(-30, 31) if m % 2 != 0])
asyms = []
for s3 in scales:
    U = d3q[f"U_{s3}"]
    d = U.shape[0] // 2
    wt = U[:d, 0] ** 2 + U[d:, 0] ** 2
    asym = np.abs(wt[odd_mu > 0].sum() - wt[odd_mu < 0].sum())
    asyms.append(asym)
ax.plot(scales, asyms, "o-", ms=3.5, color=SERIES[3])
ax.set_xlabel(r"$D_3$ scaling factor")
ax.set_ylabel("SM1 weight asymmetry")
ax.grid(True)
panel_label(ax, "(c)")

# (d) EN matrix
ax = axs[1, 1]
EN = ent["EN"]
odd = ent["odd_mu"]
im = ax.pcolormesh(odd, odd, EN, cmap="viridis", shading="auto",
                   rasterized=True)
ax.set_xlabel(r"odd mode $\mu$")
ax.set_ylabel(r"odd mode $\mu'$")
ax.set_xlim(-21, 21)
ax.set_ylim(-21, 21)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label(r"$E_N$", fontsize=7)
panel_label(ax, "(d)")

from figstyle import save
save(fig, "fig5_d3")
