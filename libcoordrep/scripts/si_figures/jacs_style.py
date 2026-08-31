"""
JACS-style figure configuration for CoordRep Supplementary Figures S2–S7.

Style: white background, no gradients, no shadows, muted palette.
Primary accent: dark orange (#C45A27).  All other series: greys / dark blue-grey.
Font: Arial/Helvetica.  Sizes per JACS guidelines.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# ── Palette ──────────────────────────────────────────
BLACK      = "#1A1A1A"
DARK_GREY  = "#4A4A4A"
MID_GREY   = "#7A7A7A"
LIGHT_GREY = "#B0B0B0"
PALE_GREY  = "#D8D8D8"
FAINT_GREY = "#ECECEC"
DARK_BLUE  = "#3A4F6F"
BLUE_GREY  = "#5C7A99"
ACCENT     = "#C45A27"       # muted / dark orange
ACCENT_LT  = "#E0915F"

# Ordered palette for multi-series (greys first, accent last)
SERIES_4 = [DARK_GREY, MID_GREY, BLUE_GREY, ACCENT]
SERIES_3 = [DARK_GREY, BLUE_GREY, ACCENT]
SERIES_2 = [DARK_GREY, ACCENT]
LAYER_COLORS = {              # L0 → L3 gradient
    "L0": DARK_GREY,
    "L1": MID_GREY,
    "L2": BLUE_GREY,
    "L3": ACCENT,
}

# ── Global RC params ─────────────────────────────────
def apply_jacs_style():
    """Apply JACS-quality matplotlib RC settings."""
    mpl.rcParams.update({
        # Font
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":          7,
        "axes.labelsize":     8,
        "axes.titlesize":     9,
        "xtick.labelsize":    7,
        "ytick.labelsize":    7,
        "legend.fontsize":    7,
        # Lines / markers
        "lines.linewidth":    1.3,
        "lines.markersize":   4,
        "scatter.edgecolors": BLACK,
        # Axes
        "axes.linewidth":     0.8,
        "axes.edgecolor":     BLACK,
        "axes.labelcolor":    BLACK,
        "xtick.color":        BLACK,
        "ytick.color":        BLACK,
        "xtick.major.width":  0.6,
        "ytick.major.width":  0.6,
        "xtick.major.size":   3,
        "ytick.major.size":   3,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        # Grid
        "axes.grid":          False,
        # Figure
        "figure.facecolor":   "white",
        "axes.facecolor":     "white",
        "savefig.facecolor":  "white",
        "savefig.dpi":        600,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.05,
        # Legend
        "legend.frameon":     False,
        "legend.borderpad":   0.3,
    })

# ── Helper: panel label ──────────────────────────────
def label_panel(ax, letter, x=-0.12, y=1.08, fontsize=10):
    """Add bold panel label (A, B, C, …) in top-left corner."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left",
            color=BLACK)

# ── Helper: light grid ───────────────────────────────
def light_grid(ax, axis="y"):
    ax.grid(axis=axis, color=FAINT_GREY, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

# ── Helper: save in three formats ────────────────────
def save_fig(fig, stem, out_dir):
    """Save figure as PDF, SVG, and 600-dpi PNG."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"{stem}.{ext}")
    print(f"  → {stem}.pdf / .svg / .png  saved to {out}")

# ── Helper: de-spine ─────────────────────────────────
def clean_axes(ax, top=False, right=False):
    ax.spines["top"].set_visible(top)
    ax.spines["right"].set_visible(right)
