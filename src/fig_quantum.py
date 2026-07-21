"""Figure 3: below-threshold quadrature physics of the bare resonator.
(a) optimal squeezing/antisqueezing spectrum (bare, critical coupling),
(b) weights of the two dominant squeezed supermodes across odd modes,
(c) squeezing versus detuning: bare output, ideal detection, and the
    3 dB critical-coupling bound.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import SERIES, OI, FULLW, panel_label

rep = np.load("q_rep.npz")
add = np.load("q_addendum.npz")
sw = np.load("q_sweep.npz", allow_pickle=True)["rows"]

fig, axs = plt.subplots(1, 3, figsize=(FULLW, 2.0))

# (a) bare spectrum
ax = axs[0]
oms = rep["oms"]
ax.plot(oms, 10 * np.log10(rep["smax_bare"]), color=SERIES[1], ls="--",
        label="antisqueezing")
ax.plot(oms, 10 * np.log10(rep["smin_bare"]), color=SERIES[0],
        label="squeezing")
ax.axhline(0, color="0.5", lw=0.6)
ax.axhline(-3.01, color="0.4", lw=0.7, ls=":")
ax.set_ylim(-4.6, 17)
ax.text(7.7, -4.3, "3 dB coupling limit", fontsize=6.2, va="bottom", ha="right")
ax.set_xlim(0, 8)
ax.set_xlabel(r"sideband frequency $\omega/\kappa$")
ax.set_ylabel("noise power (dB)")
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62),
          fontsize=6.3, borderpad=0.35)
panel_label(ax, "(a)")

# (b) supermode weights
ax = axs[1]
U = rep["sup_bare_U"]
w = rep["sup_bare_w"]
d = U.shape[0] // 2
odd_mu = np.array([m for m in range(-30, 31) if m % 2 != 0])
for k, c in zip(range(2), (SERIES[0], SERIES[1])):
    wt = U[:d, k] ** 2 + U[d:, k] ** 2
    ax.plot(odd_mu, wt, "o-", ms=2.5, lw=0.9, color=c,
            label=f"SM{k+1}: {10*np.log10(w[k]):.1f} dB")
ax.set_xlabel(r"odd mode index $\mu$")
ax.set_ylabel("supermode weight")
ax.set_xlim(-21, 21)
ax.set_ylim(-0.008, 0.202)
ax.legend(loc="upper right", fontsize=6.3, borderpad=0.35)
panel_label(ax, "(b)")

# (c) squeezing vs detuning
ax = axs[2]
z = [r["zeta0"] for r in sw]
ax.plot(z, [r["out_min_db"] for r in sw], "o-", ms=2.5, color=SERIES[0],
        label="bare output ($\\eta = 0.5$)")
ax.plot(add["zeta0"], add["ideal_db"], "s-", ms=2.5, color=SERIES[2],
        label="ideal detection ($\\eta \\to 1$)")
ax.axhline(-3.01, color="0.4", lw=0.7, ls=":")
ax.set_ylim(-26, 0)
ax.annotate("diverges at\nannihilation", xy=(10.45, -22.5), xytext=(7.3, -21.5),
            fontsize=6, ha="center",
            arrowprops=dict(arrowstyle="->", lw=0.7, color="0.3"))
ax.set_xlabel(r"pump detuning $\zeta_0$")
ax.set_ylabel("peak squeezing (dB)")
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.70),
          fontsize=6.3, borderpad=0.35)
panel_label(ax, "(c)")

from figstyle import save
save(fig, "fig3_quantum")
