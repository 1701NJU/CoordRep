"""Create the replacement Figure 2d canonicalization panel.

Outputs SVG/PDF as native vectors and a 600 dpi PNG preview. The geometry and
typography are intentionally compact enough for a one-column JACS subpanel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUTDIR = Path(__file__).resolve().parent
STEM = "Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803"


INK = "#151719"
MUTED = "#626970"
RULE = "#D6DADF"
CU = "#0D5070"
CU_EDGE = "#08364C"
PALE_BLUE = "#EAF3F7"
LEGACY = "#8064A2"
PALE_PURPLE = "#F1EDF6"
GOLD = "#A47508"
PALE_GOLD = "#FFF0C4"
GREEN = "#148653"


def add_text(ax, x, y, s, size=7.0, weight="normal", color=INK,
             ha="left", va="center", **kwargs):
    return ax.text(
        x,
        y,
        s,
        fontsize=size,
        fontfamily="Arial",
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        **kwargs,
    )


def draw_cu4_cycle(ax, cx, cy, radius=9.1):
    # The four equivalent Cu centers form the same cycle; only the red input
    # labels are arbitrary. The slight diamond rotation reads as a molecular
    # motif rather than a software graph box.
    positions = [
        (cx, cy + radius),
        (cx + radius, cy),
        (cx, cy - radius),
        (cx - radius, cy),
    ]
    for i in range(4):
        x1, y1 = positions[i]
        x2, y2 = positions[(i + 1) % 4]
        ax.plot([x1, x2], [y1, y2], color="#555B60", lw=1.35, zorder=1)

    label_offsets = [(4.0, 3.4), (4.1, -2.6), (-4.0, -4.0), (-5.0, 2.8)]
    for idx, ((x, y), (dx, dy)) in enumerate(zip(positions, label_offsets), start=1):
        ax.add_patch(Circle((x, y), 4.15, facecolor=CU, edgecolor=CU_EDGE, lw=0.8, zorder=3))
        add_text(ax, x, y - 0.1, "Cu", size=5.9, weight="bold", color="white", ha="center", zorder=4)
        ax.add_patch(Circle((x + dx, y + dy), 2.25, facecolor="white", edgecolor="#C83232", lw=0.75, zorder=5))
        add_text(ax, x + dx, y + dy, str(idx), size=5.3, weight="bold", color="#C83232", ha="center", zorder=6)


def arrow(ax, x1, y1, x2, y2, color=MUTED, lw=1.0):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=7.5,
            linewidth=lw,
            color=color,
            shrinkA=0,
            shrinkB=0,
        )
    )


def main():
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.0,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(3.45, 2.55), facecolor="white")
    ax = fig.add_axes([0.015, 0.025, 0.97, 0.955])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 74)
    ax.axis("off")

    # Panel label and concise claim.
    add_text(ax, 0.4, 71.2, "(d)", size=9.0, weight="bold")
    add_text(ax, 9.5, 71.2, "Input-order invariance for symmetric records", size=8.1, weight="bold")

    # Chemical motif and representative numbering orders.
    add_text(ax, 18.5, 63.7, r"same homogeneous Cu$_4$ cycle", size=6.7, weight="bold", ha="center")
    draw_cu4_cycle(ax, 18.5, 45.0, radius=9.0)
    add_text(ax, 18.5, 30.2, "arbitrary input numbering", size=6.0, color=MUTED, ha="center")
    add_text(ax, 18.5, 25.7, "1-2-3-4   |   3-1-4-2   |   ...", size=5.9, color=INK, ha="center")
    add_text(ax, 18.5, 21.7, "24 possible orders", size=6.1, weight="bold", color=CU, ha="center")

    # Light separator; the comparison itself is deliberately open, not a
    # programmer-style flow chart.
    ax.plot([38.2, 38.2], [20.2, 64.1], color=RULE, lw=0.85)

    # Legacy result.
    ax.add_patch(Rectangle((42.0, 60.4), 2.0, 2.0, facecolor=LEGACY, edgecolor="none"))
    add_text(ax, 46.0, 61.4, "Legacy WL-only", size=6.8, weight="bold")
    add_text(ax, 46.0, 55.0, "24", size=12.7, weight="bold", color=INK)
    add_text(ax, 54.2, 54.7, "orders", size=6.0, color=MUTED)
    arrow(ax, 64.0, 54.5, 72.2, 54.5, color=LEGACY, lw=1.25)
    for i, shade in enumerate(("#725493", "#947AB1", "#B7A5CA")):
        ax.add_patch(
            FancyBboxPatch(
                (75.2, 58.0 - i * 4.1),
                7.3,
                2.7,
                boxstyle="round,pad=0.08,rounding_size=0.45",
                facecolor=shade,
                edgecolor="none",
            )
        )
    add_text(ax, 85.4, 54.5, "3 strings", size=8.0, weight="bold", color=LEGACY)

    # Exact result.
    ax.plot([42.0, 97.0], [43.6, 43.6], color=RULE, lw=0.8)
    ax.add_patch(Rectangle((42.0, 38.7), 2.0, 2.0, facecolor=GOLD, edgecolor="none"))
    add_text(ax, 46.0, 39.7, "Exact canonical labeling", size=6.8, weight="bold")
    add_text(ax, 46.0, 32.8, "24/24", size=12.7, weight="bold", color=INK)
    add_text(ax, 60.0, 32.5, "orders", size=6.0, color=MUTED)
    arrow(ax, 68.2, 32.6, 76.3, 32.6, color=GOLD, lw=1.25)
    ax.add_patch(
        FancyBboxPatch(
            (78.4, 28.6),
            18.2,
            8.1,
            boxstyle="round,pad=0.15,rounding_size=1.0",
            facecolor=PALE_GOLD,
            edgecolor=GOLD,
            linewidth=0.9,
        )
    )
    add_text(ax, 87.5, 33.9, "1 state", size=8.3, weight="bold", color=GOLD, ha="center")
    add_text(ax, 87.5, 30.5, "CoordRep-State", size=5.2, color=GOLD, ha="center")

    # Two concise scope checks. They read like experimental totals, not code
    # diagnostics, and retain the pale-blue/gold language of the parent figure.
    ax.plot([2.0, 98.0], [16.7, 16.7], color=INK, lw=0.85)
    ax.add_patch(Rectangle((2.0, 3.0), 46.5, 10.7, facecolor=PALE_BLUE, edgecolor="none"))
    ax.add_patch(Rectangle((50.5, 3.0), 47.5, 10.7, facecolor=PALE_GOLD, edgecolor="none"))

    add_text(ax, 4.7, 11.0, "multinuclear", size=6.0, weight="bold", color=CU)
    add_text(ax, 4.7, 6.8, "55 × 100 variants", size=5.8, color=MUTED)
    add_text(ax, 46.0, 8.8, "5,500/5,500", size=7.2, weight="bold", color=GREEN, ha="right")

    add_text(ax, 53.2, 11.0, r"haptic ($\eta$)", size=6.0, weight="bold", color=GOLD)
    add_text(ax, 53.2, 6.8, "33 × 100 variants", size=5.8, color=MUTED)
    add_text(ax, 95.5, 8.8, "3,300/3,300", size=7.2, weight="bold", color=GREEN, ha="right")

    for suffix, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("png", {"dpi": 600}),
    ):
        fig.savefig(
            OUTDIR / f"{STEM}.{suffix}",
            bbox_inches="tight",
            pad_inches=0.02,
            facecolor="white",
            **kwargs,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
