"""
Shared visual style for phasewalk notebooks.

Uses seaborn's `deep` qualitative palette for categorical data
and matplotlib's perceptually uniform `cividis` for sequential data.

Usage
-----
    from _style import apply_style, COLORS, SEQ_CMAP
    apply_style()

`COLORS`  — semantic palette mapping concepts to colours.
`SEQ_CMAP` — string name of the sequential colormap for density,
             magnitude, and eigenvalue plots.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# Palette: derive named colours from seaborn's `deep` qualitative palette.
# --------------------------------------------------------------------------- #
# Colours are picked from `deep` (Bang Wong-derived, colour-blind friendly),
# but assigned to project-specific semantic roles. Keeping the mapping
# central lets us change the look in one place.
#
# `deep` palette positions (approx.):
#   0 blue     1 orange   2 green    3 red       4 purple
#   5 brown    6 pink     7 grey     8 yellow    9 cyan
DEEP = [
    "#4C72B0",  # 0 blue
    "#DD8452",  # 1 orange
    "#55A467",  # 2 green
    "#C44E52",  # 3 red
    "#8172B3",  # 4 purple
    "#937860",  # 5 brown
    "#DA8BC3",  # 6 pink
    "#8C8C8C",  # 7 grey
    "#CCB974",  # 8 yellow
    "#64B5CD",  # 9 cyan
]

COLORS = {
    # Population identity
    "excitatory": DEEP[3],   # red — fast, drives activity
    "inhibitory": DEEP[0],   # blue — slow, restrains

    # Stability classes (used in bifurcation diagrams)
    "stable":     DEEP[0],   # blue — calm, attracting
    "unstable":   DEEP[3],   # red — repelling
    "saddle":     DEEP[8],   # yellow-ochre — neither
    "limit_cycle": DEEP[3],  # red — same physical regime as unstable

    # Highlights & annotations
    "accent":     DEEP[2],   # green — for thresholds and key markers
    "muted":      "#5A6E7F",
    "grid":       "#E1E5EA",
    "background": "#FFFFFF",
}

# Sequential colormap for magnitudes, densities, eigenvalue traces.
# `cividis` is perceptually uniform AND colour-blind friendly, unlike
# the more common `viridis`.
SEQ_CMAP = "cividis"


def apply_style():
    """Configure matplotlib for publication-quality figures."""
    mpl.rcParams.update({
        # Figure
        "figure.facecolor":   COLORS["background"],
        "figure.dpi":         110,        # screen preview
        "savefig.dpi":        300,        # publication quality
        "savefig.bbox":       "tight",
        "savefig.facecolor":  COLORS["background"],
        "savefig.transparent": False,

        # Fonts — fall back gracefully across platforms
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Inter", "Source Sans Pro", "Helvetica Neue",
                               "Helvetica", "Arial", "DejaVu Sans"],
        "font.size":          11,
        "axes.titlesize":     12,
        "axes.titleweight":   "regular",
        "axes.titlepad":      14,
        "axes.labelsize":     11,
        "axes.labelpad":      6,
        "axes.labelweight":   "regular",
        "xtick.labelsize":    10,
        "ytick.labelsize":    10,
        "legend.fontsize":    10,
        "legend.title_fontsize": 10,
        "figure.titlesize":   13.5,
        "figure.titleweight": "regular",

        # Axes — minimalist, no top/right spines
        "axes.facecolor":     COLORS["background"],
        "axes.edgecolor":     COLORS["muted"],
        "axes.linewidth":     0.7,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.spines.left":   True,
        "axes.spines.bottom": True,
        "axes.grid":          True,
        "axes.axisbelow":     True,

        # Grid — soft, recedes
        "grid.color":         COLORS["grid"],
        "grid.linewidth":     0.5,
        "grid.alpha":         0.8,

        # Ticks
        "xtick.color":        COLORS["muted"],
        "ytick.color":        COLORS["muted"],
        "xtick.major.width":  0.7,
        "ytick.major.width":  0.7,
        "xtick.major.size":   3.5,
        "ytick.major.size":   3.5,
        "xtick.direction":    "out",
        "ytick.direction":    "out",

        # Legend — lighter, less boxy
        "legend.frameon":     True,
        "legend.framealpha":  0.92,
        "legend.facecolor":   COLORS["background"],
        "legend.edgecolor":   COLORS["grid"],
        "legend.borderpad":   0.5,
        "legend.handlelength": 1.6,

        # Lines
        "lines.linewidth":    1.6,
        "lines.solid_capstyle": "round",
        "lines.markersize":   5,

        # Default colormap
        "image.cmap":         SEQ_CMAP,
    })


def figure_caption(fig, text):
    """Add a small italic caption beneath the figure (publication style)."""
    fig.text(0.02, -0.015, text,
             fontsize=9, color=COLORS["muted"],
             style="italic", ha="left")
