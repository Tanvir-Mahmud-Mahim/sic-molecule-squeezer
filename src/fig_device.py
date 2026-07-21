"""Figure: device concept and FEM design (redesigned).
(a) 3D rendering of the photonic molecule, (b) mode-alignment diagram
explaining odd-mode-selective extraction, (c) cross-section |E|^2,
(d) FEM effective index, (e) integrated dispersion, (f) D2 design map.
"""
import numpy as np, json
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from mpl_toolkits.mplot3d import Axes3D  # noqa
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from figstyle import SERIES, OI, FULLW, panel_label, save
from dispersion import RingDispersion

fem = json.load(open("fem_final.json"))
sweep = json.load(open("fem_sweep.json"))
fld = np.load("mode_field.npz")

fig = plt.figure(figsize=(FULLW, 4.45), layout=None)
fig.set_constrained_layout(False)
gs = fig.add_gridspec(2, 6, height_ratios=[0.92, 1.0], hspace=0.32,
                      wspace=2.2, left=0.075, right=0.985, top=0.975,
                      bottom=0.10)

# ============ (a) 3D device rendering ============
ax = fig.add_subplot(gs[0, 0:3], projection="3d")
ax.computed_zorder = False


def torus(R, r, cx, cy, cz, nu=90, nv=24):
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, 2 * np.pi, nv)
    U, V = np.meshgrid(u, v)
    X = cx + (R + r * np.cos(V)) * np.cos(U)
    Y = cy + (R + r * np.cos(V)) * np.sin(U)
    Z = cz + 0.55 * r * np.sin(V)
    return X, Y, Z


def box(ax, x0, x1, y0, y1, z0, z1, color, alpha=1.0, zorder=1):
    V = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    faces = [[V[j] for j in f] for f in
             ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
              (2, 3, 7, 6), (1, 2, 6, 5), (0, 3, 7, 4))]
    pc = Poly3DCollection(faces, facecolor=color, edgecolor="none",
                          alpha=alpha, zorder=zorder)
    ax.add_collection3d(pc)


# substrate + buried oxide
box(ax, -9.5, 9.5, -7.2, 7.2, -2.4, -1.4, "#9aa7b1", 1.0, zorder=0)
box(ax, -9.5, 9.5, -7.2, 7.2, -1.4, -0.35, "#dbe6ee", 1.0, zorder=1)
# bus and drop waveguides
box(ax, -9.5, 9.5, 5.15, 5.85, -0.35, 0.15, "#3f4d59", 1.0, zorder=6)
box(ax, -9.5, 9.5, -5.85, -5.15, -0.35, 0.15, "#3f4d59", 1.0, zorder=2)
# main ring (R) and aux ring (R/2)
X, Y, Z = torus(3.55, 0.34, -2.5, 0.9, -0.1)
ax.plot_surface(X, Y, Z, color=OI["blue"], linewidth=0, antialiased=True,
                shade=True, alpha=1.0, zorder=5)
X, Y, Z = torus(1.78, 0.34, 3.6, -2.4, -0.1)
ax.plot_surface(X, Y, Z, color=OI["vermilion"], linewidth=0,
                antialiased=True, shade=True, alpha=1.0, zorder=5)
# pump arrow (into bus) and drop arrow
ax.quiver(-9.3, 5.5, 0.6, 3.2, 0, 0, color=OI["blue"], lw=1.6,
          arrow_length_ratio=0.16, zorder=10)
ax.quiver(5.6, -5.5, 0.6, 3.2, 0, 0, color=OI["vermilion"], lw=1.6,
          arrow_length_ratio=0.16, zorder=10)
ann = dict(fontsize=6.8, ha="center",
           bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.82",
                     alpha=0.95),
           arrowprops=dict(arrowstyle="-", lw=0.7, color="0.35",
                           shrinkA=2, shrinkB=0))
ax.text2D(0.315, 0.845, "CW pump", color=OI["blue"], fontsize=7,
          transform=ax.transAxes, ha="center")
ax.text2D(0.66, 0.025, "drop port (odd $\\mu$)", color=OI["vermilion"],
          fontsize=7, transform=ax.transAxes, ha="center")
ax.annotate("main ring $R$\n(2-FSR crystal)",
            xy=(0.285, 0.50), xytext=(0.07, 0.215),
            xycoords="axes fraction", textcoords="axes fraction", **ann)
ax.annotate("aux ring $R/2$\n(FSR $2D_1$)",
            xy=(0.745, 0.415), xytext=(0.92, 0.63),
            xycoords="axes fraction", textcoords="axes fraction", **ann)
ax.text2D(0.578, 0.47, "$J$", fontsize=7.5, ha="center",
          transform=ax.transAxes,
          bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                    alpha=0.8))
ax.set_xlim(-8, 8)
ax.set_ylim(-6.5, 6.5)
ax.set_zlim(-2.8, 3.2)
ax.set_box_aspect((1.35, 1.1, 0.42))
ax.view_init(elev=38, azim=-64)
ax.set_axis_off()
ax.text2D(0.02, 0.97, "(a)", transform=ax.transAxes, fontsize=8.5,
          fontweight="bold", va="top")

# ============ (b) mode-alignment diagram ============
ax = fig.add_subplot(gs[0, 3:6])
ax.set_xlim(-6.05, 4.75)
ax.set_ylim(-1.05, 3.75)
ax.axis("off")
mus = np.arange(-4, 5)
# main-ring resonances (top row) with mode indices
for m in mus:
    ax.plot([m, m], [2.55, 3.05], color="0.55", lw=1.2)
    ax.text(m, 3.20, f"{m:+d}" if m else "0", fontsize=6, ha="center")
# comb teeth on even modes
for m in mus[mus % 2 == 0]:
    ax.annotate("", xy=(m, 2.52), xytext=(m, 1.72),
                arrowprops=dict(arrowstyle="wedge,tail_width=0.55",
                                fc=OI["blue"], ec="none"))
# squeezed vacuum on odd modes
for m in mus[mus % 2 != 0]:
    ax.plot(m, 2.0, marker="o", ms=6, mfc="white", mec=OI["vermilion"],
            mew=1.2)
# aux-ring resonances (bottom row), aligned with odd modes only
for m in mus[mus % 2 != 0]:
    ax.plot([m, m], [0.42, 0.92], color=OI["vermilion"], lw=2.2)
    ax.annotate("", xy=(m, 1.80), xytext=(m, 0.99),
                arrowprops=dict(arrowstyle="<->", color=OI["vermilion"],
                                lw=0.9))
# kappa_P tag on the rightmost extraction arrow
ax.text(3.35, 1.38, "$\\kappa_P$", fontsize=6.8, color=OI["vermilion"],
        ha="left", va="center")
# left-hand row labels
ax.text(-6.0, 2.80, "main ring\n(FSR $D_1$)", fontsize=6.3, ha="left",
        va="center")
ax.text(-6.0, 1.40, "Purcell\nextraction", fontsize=6.3, ha="left",
        va="center", color=OI["vermilion"])
ax.text(-6.0, 0.55, "aux ring\n(FSR $2D_1$)", fontsize=6.3, ha="left",
        va="center", color=OI["vermilion"])
# legend inside, bottom band (empty region of the diagram)
h1 = ax.plot([], [], color=OI["blue"], lw=3)[0]
h2 = ax.plot([], [], "o", ms=5, mfc="white", mec=OI["vermilion"])[0]
ax.legend([h1, h2], ["comb teeth (even $\\mu$)",
                     "squeezed vacuum (odd $\\mu$)"],
          loc="upper center", fontsize=6.2, bbox_to_anchor=(0.5, 0.115),
          ncol=2, columnspacing=1.0, handletextpad=0.45,
          handlelength=1.5, frameon=False)
ax.text(0.012, 0.97, "(b)", transform=ax.transAxes, fontsize=8.5,
        fontweight="bold", va="top")

# ============ (c) cross-section mode ============
ax = fig.add_subplot(gs[1, 0:2])
E2 = fld["E2"] / fld["E2"].max()
im = ax.pcolormesh(fld["x"], fld["y"], E2, cmap="magma", rasterized=True,
                   shading="auto")
ax.add_patch(mp.Rectangle((-1.85 / 2, 0), 1.85, 0.5, fc="none", ec="w",
                          lw=0.9))
ax.annotate("4H-SiC core\n1.85 µm × 500 nm", xy=(0.7, 0.55),
            xytext=(-0.15, 1.02), color="w", fontsize=6.2, ha="center",
            arrowprops=dict(arrowstyle="-", color="w", lw=0.7))
ax.text(-1.65, -0.9, "SiO$_2$", color="w", fontsize=6.6)
ax.set_xlabel("x (µm)")
ax.set_ylabel("y (µm)")
ax.set_xlim(-1.8, 1.8)
ax.set_ylim(-1.0, 1.4)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cax = inset_axes(ax, width="4%", height="55%", loc="lower right",
                 borderpad=0.6)
cb = fig.colorbar(im, cax=cax)
cb.set_label("$|E|^2$", fontsize=6.5, labelpad=1, color="w")
cb.set_ticks([0, 1])
cb.ax.tick_params(labelsize=5.5, colors="w")
cb.outline.set_edgecolor("w")
panel_label(ax, "(c)", color="w")

# ============ (d) Dint ============
ax = fig.add_subplot(gs[1, 2:4])
d = fem["final"]
rd = RingDispersion(d["lams"], d["neff"])
R = rd.radius_for_fsr(350e9)
res = rd.resonances(R, mu_max=80)
mu = res["mu"]
ax.plot(mu, res["Dint"] / 2 / np.pi / 1e9, "o", ms=2.0, color=SERIES[0],
        label="FEM")
mm = np.linspace(-80, 80, 300)
ax.plot(mm, (res["D2"] * mm ** 2 / 2) / 2 / np.pi / 1e9, "--",
        color=SERIES[1], lw=1.1, label="$D_2\\mu^2/2$")
ax.set_xlabel("mode index $\\mu$")
ax.set_ylabel("$D_{\\rm int}/2\\pi$ (GHz)")
ax.legend(loc="upper center", fontsize=6.4, borderpad=0.35,
          handletextpad=0.5)
ax.set_ylim(-1.5, 29)
ax.grid(True)
panel_label(ax, "(d)")

# ============ (e) design map ============
ax = fig.add_subplot(gs[1, 4:6])
for j, h in enumerate([0.50, 0.60]):
    ws, d2s = [], []
    for key, dd in sweep.items():
        if abs(dd["height"] - h) < 1e-6:
            rdk = RingDispersion(dd["lams"], dd["neff"])
            Rk = rdk.radius_for_fsr(350e9)
            rk = rdk.resonances(Rk, mu_max=60)
            ws.append(dd["width"])
            d2s.append(rk["D2"] / 2 / np.pi / 1e6)
    o = np.argsort(ws)
    mk = "o-" if j == 0 else "s--"
    ax.plot(np.array(ws)[o], np.array(d2s)[o], mk, color=SERIES[j], ms=3.2,
            label=f"h = {int(h*1000)} nm")
ax.axhline(0, color="0.5", lw=0.6)
ax.plot([1.85], [res["D2"] / 2 / np.pi / 1e6], "*", ms=11,
        color=OI["vermilion"], zorder=5)
ax.annotate("selected", xy=(1.85, res["D2"] / 2 / np.pi / 1e6),
            xytext=(1.62, 2.4), fontsize=6.2, color=OI["vermilion"],
            arrowprops=dict(arrowstyle="-", lw=0.6,
                            color=OI["vermilion"]))
ax.set_xlabel("width (µm)")
ax.set_ylabel("$D_2/2\\pi$ (MHz)")
ax.legend(loc="upper right", fontsize=6.4, borderpad=0.35,
          handletextpad=0.5)
ax.set_ylim(-1.2, 27.5)
ax.grid(True)
panel_label(ax, "(e)")

save(fig, "fig1_device")
