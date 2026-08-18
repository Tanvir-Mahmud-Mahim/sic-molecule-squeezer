"""Graphical abstract: pipeline from FEM design to detected quantum noise.
Arrow positions are computed from the actual axes geometry, so the flow
arrows always sit exactly between panels.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from figstyle import SERIES, OI, FULLW, save

dat = np.load("lle_family.npz")
rep = np.load("q_rep.npz")
tw = np.load("q_tworing.npz")
N = int(dat["N"])
iz = int(np.argmin(np.abs(dat["zeta0"] - 6.5)))
psi = dat["psi"][iz]

fig = plt.figure(figsize=(FULLW, 1.52), layout=None)
fig.set_constrained_layout(False)
gs = fig.add_gridspec(1, 4, wspace=0.42, left=0.015, right=0.985,
                      top=0.775, bottom=0.245, width_ratios=[1, 1, 1, 1.12])

titles = ["FEM dispersion design", "2-FSR soliton crystal",
          "photonic-molecule extraction", "detected performance"]
axes = []

# --- 1: device top view
ax = fig.add_subplot(gs[0, 0])
axes.append(ax)
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect("equal")
ax.axis("off")
ax.add_patch(mp.Rectangle((0.3, 8.3), 9.4, 0.65, fc="#c9d7e4", ec="0.4",
                          lw=0.5))
ax.add_patch(mp.Rectangle((0.3, 1.05), 9.4, 0.65, fc="#c9d7e4", ec="0.4",
                          lw=0.5))
ax.add_patch(mp.Circle((3.6, 5.1), 2.75, fc="none", ec=OI["blue"], lw=3.4))
ax.add_patch(mp.Circle((7.95, 3.4), 1.42, fc="none", ec=OI["vermilion"],
                       lw=2.8))
ax.annotate("", xy=(2.9, 8.62), xytext=(0.7, 8.62),
            arrowprops=dict(arrowstyle="-|>", color=OI["blue"], lw=1.2))
ax.annotate("", xy=(9.55, 1.38), xytext=(8.3, 1.38),
            arrowprops=dict(arrowstyle="-|>", color=OI["vermilion"],
                            lw=1.2))
ax.text(3.6, 5.15, "4H-SiC", fontsize=5.8, ha="center", va="center")
ax.text(7.95, 3.45, "$R/2$", fontsize=5.6, ha="center", va="center")

# --- 2: crystal
ax = fig.add_subplot(gs[0, 1])
axes.append(ax)
theta = 2 * np.pi * np.arange(N) / N - np.pi
th = np.abs(np.fft.ifft(psi)) ** 2 * N
ax.plot(theta / np.pi, th, color=SERIES[0], lw=1.1)
ax.fill_between(theta / np.pi, 0, th, color=SERIES[0], alpha=0.3)
ax.set_xlim(-1, 1)
ax.set_ylim(0, 17.5)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_color("0.6")

# --- 3: squeezing spectra
ax = fig.add_subplot(gs[0, 2])
axes.append(ax)
ax.plot(rep["oms"], 10 * np.log10(rep["smin_bare"]), color="0.55", lw=1.0)
ax.plot(tw["oms"], 10 * np.log10(tw["design_smin"]), color=SERIES[2],
        lw=1.3)
ax.axhline(0, color="0.75", lw=0.5)
ax.axhline(-3.01, color="0.55", lw=0.6, ls=":")
ax.set_xlim(0, 40)
ax.set_ylim(-10.6, 1.0)
ax.set_xticks([])
ax.set_yticks([])
ax.annotate("bare (3 dB cap)", xy=(4.2, -2.6), xytext=(13.5, -4.9),
            fontsize=5.4, color="0.35",
            arrowprops=dict(arrowstyle="-", lw=0.55, color="0.5"))
ax.annotate("molecule", xy=(4.5, -7.5), xytext=(15.5, -9.7),
            fontsize=5.6, color=SERIES[2],
            arrowprops=dict(arrowstyle="-", lw=0.55, color=SERIES[2]))
for sp in ax.spines.values():
    sp.set_color("0.6")

# --- 4: badges
ax = fig.add_subplot(gs[0, 3])
axes.append(ax)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
badges = [("7.9 dB", "squeezing at 8.3 mW pump"),
          ("1.81 GHz", "squeezing bandwidth"),
          ("$E_N$ = 0.13", "entangled 30-mode lattice")]
for i, (big, small) in enumerate(badges):
    y0 = 0.695 - i * 0.345          # bottom of badge i
    ax.add_patch(mp.FancyBboxPatch((0.02, y0), 0.96, 0.30,
                                   boxstyle="round,pad=0.01",
                                   fc="#eef4fa", ec="#b9cbdc", lw=0.7))
    ax.text(0.07, y0 + 0.195, big, fontsize=7.2, fontweight="bold",
            color=OI["blue"], va="center")
    ax.text(0.07, y0 + 0.075, small, fontsize=5.2, color="0.25",
            va="center")

# --- stage titles under each panel, centered on the panel
for ax_i, t in zip(axes, titles):
    bb = ax_i.get_position()
    fig.text((bb.x0 + bb.x1) / 2, 0.10, t, fontsize=6.6, ha="center",
             color="0.15")

# --- flow arrows exactly between consecutive panels
for a, b in zip(axes[:-1], axes[1:]):
    x0 = a.get_position().x1
    x1 = b.get_position().x0
    xm = (x0 + x1) / 2
    fig.patches.append(mp.FancyArrowPatch(
        (xm - 0.013, 0.51), (xm + 0.013, 0.51),
        transform=fig.transFigure, arrowstyle="-|>",
        mutation_scale=10, color="0.35", lw=1.2))

fig.text(0.5, 0.945, "open-source pipeline:  material data "
         "$\\rightarrow$ FEM $\\rightarrow$ Lugiato-Lefever "
         "$\\rightarrow$ quantum model $\\rightarrow$ detected squeezing",
         fontsize=6.8, ha="center", color="0.25", style="italic")

save(fig, "fig0_abstract")
