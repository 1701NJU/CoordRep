from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(r"C:\Users\LUOwenlin\Desktop\PWHXXZ\coordREP\JACS")
OUT = ROOT / "rew20260817" / "CShM_Corrected_20260818"
SCAN = (
    ROOT
    / "revision_experiments"
    / "results"
    / "corrected_csd_strict_full_scan_20260723_v2_hapticityfix"
)

STEM = "Figure2_Canonicalization_Corrected_CShM"

INK = "#1E252B"
MUTED = "#5E6A73"
GRID = "#D9E0E4"
LIGHT = "#EEF2F4"
NAVY = "#2E617F"
TEAL = "#3D96A4"
ORANGE = "#D97A21"
PURPLE = "#6B55A3"
GREY = "#9DA7AD"
GOLD = "#8A6500"

EBAGAR = np.asarray([1.014856282443701, 14.312707388588121])
EBAGEV = np.asarray([0.9115136539869606, 13.581944119050016])


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.1,
            "axes.titlesize": 8.2,
            "axes.titleweight": "bold",
            "axes.labelsize": 7.2,
            "axes.linewidth": 0.65,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.fontsize": 6.2,
            "legend.frameon": False,
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )


def clean_axis(ax: mpl.axes.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.55, zorder=0)
    ax.set_axisbelow(True)


def panel_title(ax: mpl.axes.Axes, letter: str, title: str, x: float = -0.12) -> None:
    ax.text(
        x,
        1.08,
        f"({letter})",
        transform=ax.transAxes,
        fontsize=9.3,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )
    ax.set_title(title, loc="left", pad=7)


def load_corrected_ru_cn6() -> tuple[list[dict[str, str]], np.ndarray]:
    rows: list[dict[str, str]] = []
    values: list[list[float]] = []
    for path in sorted(SCAN.glob("shard_*/records_internal_licensed.csv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["metal"] != "Ru" or row["cn"] != "6":
                    continue
                refs = json.loads(row["shape_refs_json"])
                vector = json.loads(row["shape_values_json"])
                if refs != ["Oh", "TPr"] or len(vector) != 2:
                    raise RuntimeError(
                        f"Unexpected CN6 shape vector for {row['refcode']}: {refs} {vector}"
                    )
                rows.append(row)
                values.append([float(vector[0]), float(vector[1])])

    atlas = np.asarray(values, dtype=float)
    if len(rows) != 4365:
        raise RuntimeError(f"Expected 4,365 corrected Ru/CN6 rows, found {len(rows):,}")

    indexed = {row["refcode"]: atlas[i] for i, row in enumerate(rows)}
    if not np.allclose(indexed["EBAGAR"], EBAGAR, atol=1e-12):
        raise RuntimeError(f"EBAGAR mismatch: {indexed['EBAGAR']}")
    if not np.allclose(indexed["EBAGEV"], EBAGEV, atol=1e-12):
        raise RuntimeError(f"EBAGEV mismatch: {indexed['EBAGEV']}")
    return rows, atlas


def draw_panel_a(ax: mpl.axes.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(-0.05, 1.07, "(A)", fontsize=9.3, fontweight="bold", va="top")
    ax.text(0.07, 1.07, "From coordinates to one canonical record", fontsize=8.2, fontweight="bold", va="top")

    ax.text(0.02, 0.88, "Invariant geometric evidence", fontsize=7.5, fontweight="bold")
    ax.text(0.02, 0.80, "M–donor distances + donor–donor distances", fontsize=6.8, color=MUTED)

    matrix = np.asarray(
        [
            [0, 2, 3, 1, 2, 1],
            [2, 0, 1, 3, 1, 2],
            [3, 1, 0, 2, 3, 1],
            [1, 3, 2, 0, 1, 2],
            [2, 1, 3, 1, 0, 1],
            [1, 2, 1, 2, 1, 0],
        ]
    )
    colors = ["#FFFFFF", "#CFE5EB", "#73B7C5", "#3D96A4"]
    x0, y0, cell = 0.04, 0.28, 0.075
    for i in range(6):
        for j in range(6):
            ax.add_patch(
                Rectangle(
                    (x0 + j * cell, y0 + (5 - i) * cell),
                    cell * 0.92,
                    cell * 0.92,
                    facecolor=colors[int(matrix[i, j])],
                    edgecolor="white",
                    linewidth=0.5,
                )
            )
    ax.text(x0 + 3 * cell, y0 - 0.045, r"$D_{ij}$", ha="center", color=TEAL, fontsize=7.5)

    ax.add_patch(FancyArrowPatch((0.50, 0.53), (0.56, 0.53), arrowstyle="-|>", mutation_scale=12, lw=0.8, color=MUTED))
    ax.text(0.58, 0.69, "refined donor signatures", fontsize=7.2, fontweight="bold")
    ax.text(0.58, 0.625, "donor attribute: element", fontsize=5.7, color=MUTED)
    ax.text(0.58, 0.575, r"geometric $f_i$ = [$d_{M,i}$, sort($D_{i*}$)]", fontsize=6.35)
    signatures = ["N · 1.86 · 2.63 · 2.81 …", "N · 1.86 · 2.71 · 2.79 …", "O · 1.99 · 2.63 · 2.71 …"]
    for idx, label in enumerate(signatures):
        y = 0.48 - idx * 0.11
        ax.add_patch(FancyBboxPatch((0.58, y), 0.39, 0.075, boxstyle="round,pad=0.007,rounding_size=0.01", facecolor=LIGHT, edgecolor="none"))
        ax.text(0.60, y + 0.038, label, va="center", fontsize=6.4, color=MUTED)

    ax.add_patch(FancyArrowPatch((0.78, 0.23), (0.78, 0.16), arrowstyle="-|>", mutation_scale=11, lw=0.8, color=MUTED))
    ax.text(0.78, 0.105, "lexicographic order", ha="center", fontsize=7.2, fontweight="bold")
    for idx in range(6):
        x = 0.60 + idx * 0.07
        color = NAVY if idx < 3 else TEAL
        ax.add_patch(Circle((x, 0.025), 0.026, facecolor=color, edgecolor=INK, lw=0.5, clip_on=False))
        ax.text(x, 0.025, f"L{idx + 1}", ha="center", va="center", color="white", fontsize=5.6, fontweight="bold", clip_on=False)
    ax.text(0.02, 0.09, "remaining exact ties", fontsize=6.3, color=GOLD, fontweight="bold")
    ax.text(0.02, 0.02, "→ minimize the full record", fontsize=6.3, color=GOLD, fontweight="bold")


def draw_panel_b(ax: mpl.axes.Axes) -> None:
    categories = [
        "CN4 repeated\nclasses",
        "CN6 repeated\nclasses",
        "Cu$_4$ symmetric\nmetal order",
        "CN6 ligand/source\nrelabeling",
    ]
    tested = np.asarray([4, 8, 24, 100], dtype=float)
    y = np.arange(len(categories))[::-1]
    for yi, inp in zip(y, tested):
        ax.plot([1, inp], [yi, yi], color=GREY, lw=1.15, zorder=1)
        ax.scatter(inp, yi, s=31, color=GREY, edgecolor="white", linewidth=0.55, zorder=3)
        ax.scatter(1, yi, s=37, marker="D", color=TEAL, edgecolor="white", linewidth=0.55, zorder=4)
        ax.text(inp * 1.06, yi, f"{int(inp)} → 1", va="center", ha="left", fontsize=6.4)
    ax.scatter(3, y[2], s=43, facecolor="white", edgecolor=ORANGE, linewidth=1.2, zorder=5)
    ax.annotate(
        "WL-only: 3",
        xy=(3, y[2]),
        xytext=(4.25, y[2] - 0.35),
        fontsize=6.1,
        color=ORANGE,
        arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.65),
    )
    ax.text(4.25, y[2] + 0.32, "cycle/path control remained distinct", fontsize=5.9, color=MUTED)
    ax.set_xscale("log", base=2)
    ax.set_xlim(0.82, 155)
    ticks = [1, 2, 4, 8, 16, 32, 64, 128]
    ax.set_xticks(ticks, [str(v) for v in ticks])
    ax.set_yticks(y, categories)
    ax.set_xlabel("Count per audit case (log$_2$ scale)")
    panel_title(ax, "B", "Nuisance encodings converge to one exact State", x=-0.17)
    clean_axis(ax, "x")
    ax.scatter([], [], s=28, color=GREY, label="tested encodings")
    ax.scatter([], [], s=32, marker="D", color=TEAL, label="distinct exact outputs")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.25), ncol=2, handletextpad=0.4, columnspacing=1.0)


def draw_panel_c(ax: mpl.axes.Axes) -> None:
    scopes = ["Multimetal\nrecord graphs", "Haptic/π\nsite objects"]
    counts = np.asarray([5500, 3300])
    colors = [NAVY, PURPLE]
    ypos = [1, 0]
    bars = ax.barh(ypos, counts, height=0.45, color=colors, zorder=2)
    for bar, value, sub in zip(bars, counts, ["55 templates × 100", "33 templates × 100"]):
        ax.text(value - 120, bar.get_y() + bar.get_height() * 0.62, f"{value:,}/{value:,}", ha="right", va="center", color="white", fontweight="bold", fontsize=7.0)
        ax.text(value - 120, bar.get_y() + bar.get_height() * 0.27, sub, ha="right", va="center", color="white", fontsize=5.9)
    ax.set_yticks(ypos, scopes)
    ax.set_xlim(0, 6000)
    ax.set_xticks([0, 2000, 4000, 6000])
    ax.set_xlabel("Transformed record inputs")
    panel_title(ax, "C", "Expanded-scope exact recovery", x=-0.15)
    ax.text(
        0.98,
        0.985,
        "0 mismatches",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.7,
        color=TEAL,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=1.0),
    )
    clean_axis(ax, "x")


def draw_panel_d(
    atlas_ax: mpl.axes.Axes,
    relation_ax: mpl.axes.Axes,
    atlas: np.ndarray,
) -> tuple[mpl.collections.PolyCollection, int]:
    xlim = (0.0, 5.0)
    ylim = (9.0, 18.0)
    in_view = (
        (atlas[:, 0] >= xlim[0])
        & (atlas[:, 0] <= xlim[1])
        & (atlas[:, 1] >= ylim[0])
        & (atlas[:, 1] <= ylim[1])
    )
    n_view = int(in_view.sum())
    cmap = LinearSegmentedColormap.from_list("grey_density", ["#F0F3F4", "#A8B1B7", "#58646C"])
    hb = atlas_ax.hexbin(
        atlas[in_view, 0],
        atlas[in_view, 1],
        gridsize=27,
        extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
        mincnt=1,
        cmap=cmap,
        norm=LogNorm(),
        linewidths=0,
        zorder=1,
    )
    atlas_ax.plot([EBAGAR[0], EBAGEV[0]], [EBAGAR[1], EBAGEV[1]], color="white", lw=1.0, ls="--", zorder=4)
    atlas_ax.scatter(*EBAGAR, s=43, color=NAVY, edgecolor="white", linewidth=0.8, zorder=5)
    atlas_ax.scatter(*EBAGEV, s=43, color=ORANGE, edgecolor="white", linewidth=0.8, zorder=5)
    atlas_ax.annotate("fac\n1.015, 14.313", EBAGAR, xytext=(9, 7), textcoords="offset points", fontsize=5.9, color=NAVY, fontweight="bold", linespacing=1.0)
    atlas_ax.annotate("mer\n0.912, 13.582", EBAGEV, xytext=(9, -15), textcoords="offset points", fontsize=5.9, color=ORANGE, fontweight="bold", linespacing=1.0)
    atlas_ax.set_xlim(*xlim)
    atlas_ax.set_ylim(*ylim)
    atlas_ax.set_xticks([0, 1, 2, 3, 4, 5])
    atlas_ax.set_yticks([9, 12, 15, 18])
    atlas_ax.set_xlabel(r"$S$(Oh)")
    atlas_ax.set_ylabel(r"$S$(TPr)")
    panel_title(atlas_ax, "D", "Corrected CShM atlas", x=-0.23)
    atlas_ax.text(
        0.02,
        0.98,
        f"strict Ru/CN = 6: {len(atlas):,} total\n{n_view:,} ({100*n_view/len(atlas):.1f}%) in view",
        transform=atlas_ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.6,
        color=MUTED,
        linespacing=1.05,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.86, pad=1.4),
    )
    clean_axis(atlas_ax)

    pair_classes = [r"N$_L$–Cl", r"N$_L$–N$_{NO}$", r"N$_L$–N$_L$", r"N$_{NO}$–Cl"]
    fac_counts = np.asarray([2, 1, 0, 0])
    mer_counts = np.asarray([1, 0, 1, 1])
    yy = np.arange(len(pair_classes))[::-1]
    height = 0.32
    relation_ax.barh(yy + height / 2, fac_counts, height, color=NAVY, label="EBAGAR fac", zorder=2)
    relation_ax.barh(yy - height / 2, mer_counts, height, color=ORANGE, label="EBAGEV mer", zorder=2)
    for y, fa, me in zip(yy, fac_counts, mer_counts):
        if fa:
            relation_ax.text(fa + 0.06, y + height / 2, str(int(fa)), va="center", fontsize=5.6, color=NAVY)
        if me:
            relation_ax.text(me + 0.06, y - height / 2, str(int(me)), va="center", fontsize=5.6, color=ORANGE)
    relation_ax.set_yticks(yy, pair_classes)
    relation_ax.set_xlim(0, 2.45)
    relation_ax.set_ylim(-0.55, 4.25)
    relation_ax.set_xticks([0, 1, 2])
    relation_ax.set_xlabel("trans pairs")
    relation_ax.set_title("Trans relations", loc="left", pad=7)
    relation_ax.tick_params(axis="y", labelsize=5.8, pad=2)
    relation_ax.text(0.02, 0.99, "■ EBAGAR fac", transform=relation_ax.transAxes, ha="left", va="top", fontsize=5.4, color=NAVY)
    relation_ax.text(0.02, 0.92, "■ EBAGEV mer", transform=relation_ax.transAxes, ha="left", va="top", fontsize=5.4, color=ORANGE)
    clean_axis(relation_ax, "x")
    return hb, n_view


def write_source_data(
    hb: mpl.collections.PolyCollection,
    n_atlas: int,
    n_view: int,
) -> None:
    source_path = OUT / f"{STEM}_source_data.csv"
    with source_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "series", "category", "value", "unit", "source_note"])
        for category, value in zip(
            ["CN4 repeated classes", "CN6 repeated classes", "Cu4 symmetric metal order", "CN6 global relabeling"],
            [4, 8, 24, 100],
        ):
            writer.writerow(["B", "tested encodings", category, value, "inputs", "exact canonicalization audit"])
            writer.writerow(["B", "distinct exact outputs", category, 1, "State", "exact canonicalization audit"])
        writer.writerow(["B", "preliminary WL-only", "Cu4 symmetric metal order", 3, "serializations", "24 input permutations"])
        writer.writerow(["C", "exact recovery", "multimetal record graphs", 5500, "of 5500", "55 curated templates x 100"])
        writer.writerow(["C", "exact recovery", "haptic/pi site objects", 3300, "of 3300", "33 curated templates x 100"])
        writer.writerow(["D", "corrected strict atlas", "Ru-centered CN=6 total", n_atlas, "records", "four-shard corrected CSD strict scan"])
        writer.writerow(["D", "corrected strict atlas", "records in plotted window", n_view, "records", "0<=S(Oh)<=5 and 9<=S(TPr)<=18"])
        writer.writerow(["D", "corrected strict atlas", "records outside plotted window", n_atlas - n_view, "records", "not displayed; included in total"])
        writer.writerow(["D", "CShM", "EBAGAR fac S(Oh)", EBAGAR[0], "CShM", "corrected exact coordinate-space calculation"])
        writer.writerow(["D", "CShM", "EBAGAR fac S(TPr)", EBAGAR[1], "CShM", "corrected exact coordinate-space calculation"])
        writer.writerow(["D", "CShM", "EBAGEV mer S(Oh)", EBAGEV[0], "CShM", "corrected exact coordinate-space calculation"])
        writer.writerow(["D", "CShM", "EBAGEV mer S(TPr)", EBAGEV[1], "CShM", "corrected exact coordinate-space calculation"])
        relation_rows = [
            ("EBAGAR fac", "N_L-Cl", 2),
            ("EBAGAR fac", "N_L-N_NO", 1),
            ("EBAGAR fac", "N_L-N_L", 0),
            ("EBAGAR fac", "N_NO-Cl", 0),
            ("EBAGEV mer", "N_L-Cl", 1),
            ("EBAGEV mer", "N_L-N_NO", 0),
            ("EBAGEV mer", "N_L-N_L", 1),
            ("EBAGEV mer", "N_NO-Cl", 1),
        ]
        for series, category, value in relation_rows:
            writer.writerow(["D", series, category, value, "trans pairs", "CoordRep relation field"])

    hex_path = OUT / f"{STEM}_Ru_CN6_hexbin.csv"
    offsets = hb.get_offsets()
    counts = hb.get_array()
    with hex_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["S_Oh_bin_center", "S_TPr_bin_center", "count"])
        for (x, y), count in zip(offsets, counts):
            writer.writerow([f"{float(x):.8f}", f"{float(y):.8f}", int(round(float(count)))])


def write_caption_and_readme(n_atlas: int, n_view: int) -> None:
    caption = (
        "Figure 2. Exact canonicalization removes input conventions while retaining coordination stereochemistry. "
        "(A) Rotation-invariant metal–donor and donor–donor distances define refined donor signatures; lexicographic refinement fixes donor order, and residual exact ties are resolved by minimizing the full record. "
        "(B) Exact tie resolution converged to one CoordRep-State for all four CN = 4 repeated-class labelings, all eight CN = 6 repeated-class labelings, the 24 orderings of a symmetric Cu4 cycle, and 100/100 random global ligand-order and source-relabeling trials; the cycle remained distinct from the corresponding path graph. "
        "(C) Canonicalization reproduced one record for every transformed input from 55 curated multimetal templates (5,500/5,500) and 33 curated haptic templates (3,300/3,300). "
        f"(D) Corrected coordinate-space CShM values were recomputed for a strict mononuclear, atom-resolved η¹ cohort of {n_atlas:,} Ru-centered CN = 6 CSD records. "
        f"The plotted window contains {n_view:,}/{n_atlas:,} records ({100*n_view/n_atlas:.1f}%). "
        "EBAGAR (fac) has [S(Oh), S(TPr)] = [1.014856, 14.312707], whereas EBAGEV (mer) has [0.911514, 13.581944]. "
        "Their explicit trans-donor inventories remain distinct despite a shared connectivity-level identity."
    )
    (OUT / f"{STEM}_caption.txt").write_text(caption + "\n", encoding="utf-8")

    readme = f"""# Corrected Figure 2 evidence package

This package replaces the old Figure 2D CShM atlas and the values formerly reported for EBAGAR/EBAGEV.

## What changed

- The old 4,449-record Ru/CN=6 atlas was retired because it used the superseded CShM objective.
- The background was rebuilt from the four completed shards in `corrected_csd_strict_full_scan_20260723_v2_hapticityfix`.
- The corrected strict cohort contains **{n_atlas:,}** Ru-centered CN=6 records; **{n_view:,} ({100*n_view/n_atlas:.1f}%)** fall in the plotted window.
- EBAGAR: `S(Oh)=1.0148562824`, `S(TPr)=14.3127073886`.
- EBAGEV: `S(Oh)=0.9115136540`, `S(TPr)=13.5819441191`.
- The fac/mer claim is supported by the explicit trans-donor relation inventory, not by a large separation in CShM space.

## Scope and provenance

The atlas is the strict mononuclear, nondisordered, nonpolymeric, atom-resolved eta1 core emitted by CoordRep 1.1.2rc2 from the April 2025 CSD scan. CShM values use the corrected exact coordinate-space calculation. The row-level CSD-derived files are licensed internal derivatives and are not copied into this package. `*_Ru_CN6_hexbin.csv` contains only aggregate bin counts.

## Rebuild

Run `py -3.12 generate_figure2_corrected_cshm.py` from this directory (or pass its full path). The script checks the expected cohort size and both highlighted CShM vectors before writing outputs.

## Outputs

- `{STEM}.svg`: editable vector artwork.
- `{STEM}.pdf`: one-page vector PDF.
- `{STEM}_600dpi.png`: 600 dpi raster for Word/layout proofing.
- `{STEM}_source_data.csv`: compact panel data and exact highlighted values.
- `{STEM}_Ru_CN6_hexbin.csv`: aggregate corrected atlas bins (no row-level CSD records).
- `{STEM}_caption.txt`: replacement caption.
"""
    (OUT / f"{STEM}_README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    configure_matplotlib()
    OUT.mkdir(parents=True, exist_ok=True)
    _rows, atlas = load_corrected_ru_cn6()

    fig = plt.figure(figsize=(178 / 25.4, 118 / 25.4))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.02, 0.98],
        width_ratios=[0.96, 1.18],
        left=0.095,
        right=0.975,
        top=0.95,
        bottom=0.105,
        wspace=0.37,
        hspace=0.55,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    draw_panel_a(ax_a)

    ax_b = fig.add_subplot(outer[0, 1])
    draw_panel_b(ax_b)

    ax_c = fig.add_subplot(outer[1, 0])
    draw_panel_c(ax_c)

    dgrid = outer[1, 1].subgridspec(1, 2, width_ratios=[1.18, 0.82], wspace=0.48)
    ax_d1 = fig.add_subplot(dgrid[0, 0])
    ax_d2 = fig.add_subplot(dgrid[0, 1])
    hb, n_view = draw_panel_d(ax_d1, ax_d2, atlas)

    svg = OUT / f"{STEM}.svg"
    pdf = OUT / f"{STEM}.pdf"
    png = OUT / f"{STEM}_600dpi.png"
    fig.savefig(svg)
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    write_source_data(hb, len(atlas), n_view)
    write_caption_and_readme(len(atlas), n_view)
    plt.close(fig)

    # Short aliases make the corrected asset easy to identify beside the prior
    # Illustrator-exported AFIG2.svg without changing the original file.
    shutil.copyfile(svg, OUT / "AFIG2_corrected.svg")
    shutil.copyfile(pdf, OUT / "AFIG2_corrected.pdf")
    shutil.copyfile(png, OUT / "AFIG2_corrected_600dpi.png")

    print(f"wrote {svg}")
    print(f"atlas records: {len(atlas):,}; plotted: {n_view:,}")


if __name__ == "__main__":
    main()
