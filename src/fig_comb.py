"""Figure 2: the 2-FSR soliton-crystal mean field.
(a) intracavity intensity profile, (b) comb spectrum (even/odd),
(c) existence diagram: peak intensity vs detuning with stationary and
    breathing windows, (d) D3-induced repetition-rate shift vs detuning.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import SERIES, OI, FULLW, panel_label

dat = np.load("lle_family.npz")
N = int(dat["N"])
kappa = float(dat["kappa"])
mu_fft = np.fft.fftfreq(N, d=1.0 / N).astype(int)
zeta0 = dat["zeta0"]
resid = dat["resid"]
crystal = dat["crystal"]
vs = dat["v"]
psis = dat["psi"]

iz = int(np.argmin(np.abs(zeta0 - 6.5)))
psi = psis[iz]

fig, axs = plt.subplots(2, 2, figsize=(FULLW, 3.55))

# (a) intensity profile
ax = axs[0, 0]
theta = 2 * np.pi * np.arange(N) / N - np.pi
th = np.abs(np.fft.ifft(psi)) ** 2 * N
ax.plot(theta / np.pi, np.fft.fftshift(th) if False else th, color=SERIES[0])
ax.fill_between(theta / np.pi, 0, th, color=SERIES[0], alpha=0.25)
ax.set_xlabel(r"azimuthal angle $\theta/\pi$")
ax.set_ylabel(r"$|\psi(\theta)|^2$")
ax.set_xlim(-1, 1)
ax.text(0.035, 0.80, rf"$\zeta_0 = {zeta0[iz]:.1f}$", transform=ax.transAxes,
        fontsize=7.5, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.85"))
panel_label(ax, "(a)")

# (b) comb spectrum
ax = axs[0, 1]
p = np.abs(psi) ** 2 / N
pdb = 10 * np.log10(p / p[0] + 1e-30)
mus = np.arange(-40, 41)
vals = [pdb[np.where(mu_fft == m)[0][0]] for m in mus]
even = np.array([m % 2 == 0 for m in mus])
ax.vlines(mus[even], -245, np.array(vals)[even], color=SERIES[0], lw=1.4,
          label="even $\\mu$ (comb)")
ax.plot(mus[~even], np.array(vals)[~even], "v", ms=2.5, color=SERIES[3],
        label="odd $\\mu$ (numerical floor)")
ax.set_ylim(-245, 8)
ax.set_xlabel("mode index $\\mu$")
ax.set_ylabel("power (dB)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
          fontsize=6.4, frameon=False, borderpad=0.1,
          columnspacing=1.0, handletextpad=0.5)
panel_label(ax, "(b)")

# (c) existence diagram
ax = axs[1, 0]
peak = [np.max(np.abs(np.fft.ifft(ps)) ** 2 * N) for ps in psis]
stat = (resid < 1e-2) & crystal
brth = (resid >= 1e-2) & crystal
ax.plot(zeta0[stat], np.array(peak)[stat], "o", ms=3, color=SERIES[0],
        label="stationary crystal")
ax.plot(zeta0[brth], np.array(peak)[brth], "s", ms=3, mfc="none",
        color=SERIES[1], label="breathing")
ax.plot(zeta0[~crystal], np.array(peak)[~crystal], "x", ms=4,
        color="0.6", label="lost")
zmax = np.pi ** 2 * float(dat["f"]) ** 2 / 8
ax.axvline(zmax, color="0.5", lw=0.7, ls=":")
ax.annotate(r"$\pi^2 f^2/8$", xy=(zmax, 21.3), xytext=(9.35, 22.6),
            fontsize=6.5, ha="center",
            arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
ax.set_ylim(4, 24.5)
ax.set_xlabel(r"pump detuning $\zeta_0 = 2\delta_0/\kappa$")
ax.set_ylabel(r"peak $|\psi|^2$")
ax.legend(loc="lower right", fontsize=6.4, borderpad=0.3)
panel_label(ax, "(c)")

# (d) rep-rate shift
ax = axs[1, 1]
vhz = vs * (kappa / 2) / (2 * np.pi) / 1e6   # v in normalized units -> MHz
ax.plot(zeta0[stat], vhz[stat], "o-", ms=3, color=SERIES[2])
ax.set_xlabel(r"pump detuning $\zeta_0$")
ax.set_ylabel(r"$\delta f_{\rm rep}$ (MHz)")
ax.grid(True)
panel_label(ax, "(d)")

from figstyle import save
save(fig, "fig2_comb")
