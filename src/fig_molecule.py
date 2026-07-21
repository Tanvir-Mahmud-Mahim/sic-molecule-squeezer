"""Figure 4: photonic-molecule outcoupling of the odd-mode lattice.
(a) drop-port squeezing spectra for increasing Purcell rate,
(b) peak squeezing vs kappa_P (with escape efficiency annotated) at
    fixed pump; (c) squeezing bandwidth vs kappa_P; (d) squeezing vs
    detuning for kP = 10 and 20 (detuning as a squeezing reservoir).
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import SERIES, OI, FULLW, panel_label

rep = np.load("q_rep.npz")
add = np.load("q_addendum.npz")

fig, axs = plt.subplots(2, 2, figsize=(FULLW, 3.6))

# (a) spectra
ax = axs[0, 0]
omsm = rep["omsm"]
sel = [0.5, 5.0, 10.0, 20.0, 50.0]
for j, kp in enumerate(sel):
    lab = "bare (0.5)" if kp == 0.5 else f"{kp:g}"
    ax.plot(omsm, 10 * np.log10(rep[f"smin_{kp}"]), color=SERIES[j % 6],
            lw=1.0, label=lab)
ax.axhline(0, color="0.5", lw=0.6)
ax.set_xlim(0, 45)
ax.set_xlabel(r"sideband frequency $\omega/\kappa$")
ax.set_ylabel("squeezing (dB)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=5,
          fontsize=5.9, title=r"$\kappa_P/\kappa$:", title_fontsize=6.2,
          frameon=False, borderpad=0.1, columnspacing=0.7,
          handletextpad=0.35, handlelength=1.2)
panel_label(ax, "(a)")

# (b) peak squeezing vs kP
ax = axs[0, 1]
kps = rep["kps"]
ax.semilogx(kps, rep["mol_min_db"], "o-", ms=3.5, color=SERIES[0])
for kp, sv, e, off in zip(kps, rep["mol_min_db"], rep["mol_eta"],
                          [None, None, (8, 4), None, (0, 9), None, None,
                           (-14, -13)]):
    if off is not None:
        ax.annotate(f"$\\eta$ = {e:.2f}", (kp, sv),
                    textcoords="offset points", xytext=off, fontsize=6.4,
                    ha="center")
ax.set_xlabel(r"Purcell rate $\kappa_P/\kappa$")
ax.set_ylabel("peak squeezing (dB)")
ax.grid(True, which="both")
panel_label(ax, "(b)")

# (c) bandwidth
ax = axs[1, 0]
ax.semilogx(kps, rep["mol_bw"] * 0.1289, "s-", ms=3.5, color=SERIES[2])
ax.set_xlabel(r"Purcell rate $\kappa_P/\kappa$")
ax.set_ylabel("squeezing bandwidth (GHz)")
ax.grid(True, which="both")
panel_label(ax, "(c)")

# (d) squeezing vs detuning for molecule
ax = axs[1, 1]
ax.plot(add["zeta0"], add["mol10_db"], "o-", ms=2.5, color=SERIES[0],
        label=r"$\kappa_P = 10\kappa$")
ax.plot(add["zeta0"], add["mol20_db"], "s-", ms=2.5, color=SERIES[1],
        label=r"$\kappa_P = 20\kappa$")
ax.plot(add["zeta0"], add["ideal_db"], ":", color="0.4",
        label=r"ideal ($\eta \to 1$)")
ax.set_ylim(-26, 0)
ax.set_xlabel(r"pump detuning $\zeta_0$")
ax.set_ylabel("peak squeezing (dB)")
ax.legend(fontsize=6.2, loc="upper right", borderpad=0.3)
panel_label(ax, "(d)")

from figstyle import save
save(fig, "fig4_molecule")
