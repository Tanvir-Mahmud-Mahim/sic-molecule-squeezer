"""Shared publication figure style (Optica universal template, 5.9 in text).

Categorical colors: Okabe-Ito (CVD-safe), fixed assignment order.
Sequential maps: single-hue perceptual (viridis/magma family).
All figures use constrained_layout to guarantee no overlaps.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OI = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermilion": "#D55E00",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "black": "#000000",
}
SERIES = [OI["blue"], OI["orange"], OI["green"], OI["vermilion"],
          OI["sky"], OI["purple"]]

plt.rcParams.update({
    "font.size": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 0.7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "legend.edgecolor": "0.85",
    "legend.fancybox": False,
    "lines.linewidth": 1.3,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "grid.alpha": 0.22,
    "grid.linewidth": 0.4,
    "pdf.fonttype": 42,
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.06,
    "figure.constrained_layout.w_pad": 0.06,
})

FULLW = 5.9   # inches, full text width
COLW = 3.25


def panel_label(ax, s, loc="in", color="k"):
    """Panel label inside the axes, top-left, with a soft white backing."""
    if loc == "in":
        ax.text(0.03, 0.965, s, transform=ax.transAxes, fontsize=8.5,
                fontweight="bold", va="top", ha="left", color=color,
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="none", alpha=0.75), zorder=20)
    else:
        ax.text(-0.14, 1.02, s, transform=ax.transAxes, fontsize=8.5,
                fontweight="bold", va="bottom", ha="left", color=color,
                zorder=20)


def save(fig, name):
    # PDF without bbox-tight: tight cropping misplaces rasterized artists
    # in some PDF viewers; margins are handled by layout instead.
    fig.savefig(f"../figures/{name}.pdf")
    fig.savefig(f"../figures/{name}.png")
    print("saved", name)
