"""Figure 6: explicit two-ring model of the photonic molecule.
(a) heater-aligned resonance mismatch across the odd lattice for three
    auxiliary-radius errors,
(b) resulting mode-dependent Purcell rate at the design point,
(c) drop-port squeezing spectra: adiabatic model versus the full two-ring
    model for increasing auxiliary linewidth,
(d) peak squeezing versus auxiliary linewidth (adiabatic-elimination
    convergence) with radius-error tolerance at the design point.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import SERIES, OI, FULLW, panel_label, save

aux = np.load("aux_ring.npz")
tw = np.load("q_tworing.npz")
kappa = float(tw["kappa"])
mu = np.asarray(aux["mu_odd"])
kaux8 = 4 * float(aux["kaux"])                    # rad/s, 8 GHz design
KP0 = 10.0

fig, axs = plt.subplots(2, 2, figsize=(FULLW, 3.05))

# (a) mismatch profiles
ax = axs[0, 0]
labels = {"0nm": "ideal radius", "5nm": r"$\delta R = 5$ nm",
          "20nm": r"$\delta R = 20$ nm"}
for k, (name, lab) in enumerate(labels.items()):
    ax.plot(mu, np.asarray(aux[f"dm_{name}"]) / 2 / np.pi / 1e9, "o-",
            ms=2.6, lw=1.0, color=SERIES[k], label=lab)
ax.axhspan(-kaux8 / 2 / 2 / np.pi / 1e9, kaux8 / 2 / 2 / np.pi / 1e9,
           color="0.85", alpha=0.45, lw=0, zorder=0)
ax.text(28.8, -3.35, r"$\pm\kappa_{\rm aux}/2$", fontsize=6.2,
        color="0.35", ha="right", va="center")
ax.set_xlabel(r"odd mode index $\mu$")
ax.set_ylabel(r"mismatch $\Delta_\mu/2\pi$ (GHz)")
ax.set_xlim(-31, 31)
ax.legend(fontsize=6.2, loc="upper left", bbox_to_anchor=(0.015, 0.86),
          borderpad=0.3, handletextpad=0.5)
ax.grid(True)
panel_label(ax, "(a)")

# (b) kappa_P(mu) at the design point
ax = axs[0, 1]
for k, (name, lab) in enumerate(labels.items()):
    dm = np.asarray(aux[f"dm_{name}"])
    kpmu = KP0 / (1 + (2 * dm / kaux8) ** 2)
    ax.plot(mu, kpmu, "o-", ms=2.6, lw=1.0, color=SERIES[k], label=lab)
ax.set_xlabel(r"odd mode index $\mu$")
ax.set_ylabel(r"$\kappa_P(\mu)/\kappa$")
ax.set_xlim(-31, 31)
ax.set_ylim(0, 11.2)
ax.legend(fontsize=6.2, loc="lower center", borderpad=0.3,
          handletextpad=0.5)
ax.grid(True)
panel_label(ax, "(b)")

# (c) spectra: adiabatic vs full for growing kappa_aux
ax = axs[1, 0]
oms = tw["oms"]
ax.plot(oms, 10 * np.log10(tw["smin_ad_10"]), color="0.2", lw=1.4,
        ls="--", label="adiabatic")
ax.plot(oms, 10 * np.log10(tw["smin_full2_10"]), color=SERIES[4], lw=1.1,
        label=r"full, $\kappa_{\rm aux}/2\pi = 2$ GHz")
ax.plot(oms, 10 * np.log10(tw["smin_full8_10"]), color=SERIES[0], lw=1.3,
        label=r"full, $\kappa_{\rm aux}/2\pi = 8$ GHz")
ax.plot(oms, 10 * np.log10(tw["design_smin_20nm"]), color=SERIES[3],
        lw=1.1, label=r"full, $\delta R = 20$ nm")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_xlim(0, 40)
ax.set_ylim(-9.4, 0.8)
ax.set_xlabel(r"sideband frequency $\omega/\kappa$")
ax.set_ylabel("squeezing (dB)")
ax.legend(fontsize=6.2, loc="lower right", borderpad=0.3)
ax.grid(True)
panel_label(ax, "(c)")

# (d) convergence + tolerance
ax = axs[1, 1]
ck = np.asarray(tw["conv_kaux"]) * kappa / 2 / np.pi / 1e9
cd = np.asarray(tw["conv_db"])
ax.semilogx(ck, -cd, "o-", ms=3.5, color=SERIES[2],
            label="full two-ring model")
ax.axhline(-float(tw["peak_ad_10"]), color="0.2", lw=1.0, ls="--",
           label="adiabatic limit")
for name, c, lab in (("5nm", SERIES[1], r"$\delta R = 5$ nm"),
                     ("20nm", SERIES[3], r"$\delta R = 20$ nm")):
    ax.plot([8], [-float(tw[f"design_real_{name}"])], "s", ms=4.5,
            color=c, label=lab)
ax.set_xlabel(r"auxiliary linewidth $\kappa_{\rm aux}/2\pi$ (GHz)")
ax.set_ylabel("peak squeezing (dB)")
ax.set_ylim(6.4, 8.8)
ax.legend(fontsize=6.2, loc="lower right", borderpad=0.3)
ax.grid(True, which="both")
panel_label(ax, "(d)")

save(fig, "fig6_tworing")
