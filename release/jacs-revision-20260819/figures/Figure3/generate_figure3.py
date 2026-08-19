from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from PIL import Image


ROOT = Path(r"C:\Users\LUOwenlin\Desktop\PWHXXZ\coordREP\JACS")
SOURCE = (
    ROOT
    / "revision_experiments"
    / "results"
    / "corrected_cshm_xtb_atlas_v1"
    / "corrected_shape_atlas.csv"
)
SOURCE_SUMMARY = SOURCE.with_name("SUMMARY.json")
TRANSITIONS = SOURCE.with_name("frozen_to_corrected_transitions.csv")
OUT = ROOT / "rew20260817" / "CShM_Corrected_20260818"

STEM = "Figure3_Corrected_CShM_Atlas_JACS_20260818"

INK = "#1E252B"
MUTED = "#5C6870"
GRID = "#D9E0E4"
NAVY = "#2E617F"
TEAL = "#2F8792"
ORANGE = "#D97A21"
PURPLE = "#6B55A3"
RED = "#B84B4B"

# A restrained, print-friendly sequential map with a white low-count end.
DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "coordrep_density", ["#EEF2F4", "#A8CBD4", "#2F8792", "#1E4D66"]
)

PAIR_SPECS = {
    4: {
        "reference_a": "SP",
        "reference_b": "Td",
        "column_a": "cshm_SP",
        "column_b": "cshm_Td",
        "color_a": NAVY,
        "color_b": ORANGE,
    },
    5: {
        "reference_a": "SPY",
        "reference_b": "TBP",
        "column_a": "cshm_SPY",
        "column_b": "cshm_TBP",
        "color_a": TEAL,
        "color_b": PURPLE,
    },
    6: {
        "reference_a": "Oh",
        "reference_b": "TPr",
        "column_a": "cshm_Oh",
        "column_b": "cshm_TPr",
        "color_a": NAVY,
        "color_b": RED,
    },
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.3,
            "axes.titlesize": 8.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 7.1,
            "axes.linewidth": 0.65,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.8,
            "ytick.major.size": 2.8,
            "figure.dpi": 150,
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
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.55, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax: mpl.axes.Axes, letter: str) -> None:
    ax.text(
        -0.17,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=9.2,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_source(frame: pd.DataFrame, summary: dict[str, object]) -> None:
    required = {
        "mol_id",
        "cn",
        "best_shape",
        "best_cshm",
        "second_shape",
        "shape_margin",
        "boundary_margin_lt_1",
        "cshm_SP",
        "cshm_Td",
        "cshm_TBP",
        "cshm_SPY",
        "cshm_Oh",
        "cshm_TPr",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing source columns: {missing}")
    if int(summary["input_records"]) != 48_057:
        raise ValueError("Unexpected cohort size; Figure 3 is locked to 48,057 records")
    if int(summary["shape_supported_records"]) != 35_485:
        raise ValueError("Unexpected supported-shape denominator")
    if len(frame) != 35_485:
        raise ValueError(f"Atlas row count is {len(frame):,}, expected 35,485")


def write_panel_source(frame: pd.DataFrame) -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    for cn, spec in PAIR_SPECS.items():
        subset = frame.loc[frame["cn"] == cn].copy()
        records.append(
            pd.DataFrame(
                {
                    "mol_id": subset["mol_id"],
                    "cn": cn,
                    "reference_a": spec["reference_a"],
                    "reference_b": spec["reference_b"],
                    "S_reference_a": subset[spec["column_a"]],
                    "S_reference_b": subset[spec["column_b"]],
                    "best_shape": subset["best_shape"],
                    "best_cshm": subset["best_cshm"],
                    "second_shape": subset["second_shape"],
                    "shape_margin": subset["shape_margin"],
                    "boundary_margin_lt_1": subset["boundary_margin_lt_1"],
                }
            )
        )
    panel_data = pd.concat(records, ignore_index=True)
    panel_data.to_csv(OUT / "Figure3_panel_source.csv", index=False)
    return panel_data


def build_summary(panel_data: pd.DataFrame) -> pd.DataFrame:
    transition = pd.read_csv(TRANSITIONS)
    rows: list[dict[str, object]] = []
    for cn, spec in PAIR_SPECS.items():
        subset = panel_data.loc[panel_data["cn"] == cn]
        transitions_cn = transition.loc[transition["cn"] == cn]
        retained = int(
            transitions_cn.loc[
                transitions_cn["frozen_best_shape"]
                == transitions_cn["corrected_best_shape"],
                "n",
            ].sum()
        )
        n = len(subset)
        boundary_n = int(subset["boundary_margin_lt_1"].sum())
        rows.append(
            {
                "scope": f"CN{cn}",
                "n": n,
                "reference_a": spec["reference_a"],
                "nearest_reference_a_n": int(
                    (subset["best_shape"] == spec["reference_a"]).sum()
                ),
                "reference_b": spec["reference_b"],
                "nearest_reference_b_n": int(
                    (subset["best_shape"] == spec["reference_b"]).sum()
                ),
                "boundary_margin_lt_1_n": boundary_n,
                "boundary_margin_lt_1_percent": 100.0 * boundary_n / n,
                "internal_frozen_label_retained_n": retained,
                "internal_frozen_label_retained_percent": 100.0 * retained / n,
            }
        )
    cn456 = panel_data.loc[panel_data["cn"].isin(PAIR_SPECS)]
    boundary_n = int(cn456["boundary_margin_lt_1"].sum())
    rows.append(
        {
            "scope": "CN4-CN6",
            "n": len(cn456),
            "reference_a": "",
            "nearest_reference_a_n": "",
            "reference_b": "",
            "nearest_reference_b_n": "",
            "boundary_margin_lt_1_n": boundary_n,
            "boundary_margin_lt_1_percent": 100.0 * boundary_n / len(cn456),
            "internal_frozen_label_retained_n": int(
                sum(row["internal_frozen_label_retained_n"] for row in rows)
            ),
            "internal_frozen_label_retained_percent": 100.0
            * sum(row["internal_frozen_label_retained_n"] for row in rows)
            / len(cn456),
        }
    )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "Figure3_summary_and_internal_QC.csv", index=False)
    return result


def draw_shape_panel(
    ax: mpl.axes.Axes,
    panel_data: pd.DataFrame,
    cn: int,
    letter: str,
) -> mpl.collections.PolyCollection:
    spec = PAIR_SPECS[cn]
    subset = panel_data.loc[panel_data["cn"] == cn]
    x = np.log10(subset["S_reference_a"].to_numpy(dtype=float) + 0.1)
    y = np.log10(subset["S_reference_b"].to_numpy(dtype=float) + 0.1)
    extent = (-1.05, 1.72, -1.05, 1.72)
    density = ax.hexbin(
        x,
        y,
        gridsize=36,
        extent=extent,
        mincnt=1,
        cmap=DENSITY_CMAP,
        norm=LogNorm(vmin=1, vmax=1_100),
        linewidths=0,
        zorder=2,
    )
    if int(np.asarray(density.get_array()).sum()) != len(subset):
        raise ValueError(f"CN{cn} hexbin extent omitted one or more records")
    ax.plot(
        [extent[0], extent[1]],
        [extent[0], extent[1]],
        color=MUTED,
        linewidth=0.75,
        linestyle=(0, (3.0, 2.2)),
        zorder=3,
    )
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.set_aspect("equal", adjustable="box")
    ticks = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.grid(color=GRID, linewidth=0.42, zorder=0)
    ax.set_axisbelow(True)
    clean_axis(ax)
    ref_a = str(spec["reference_a"])
    ref_b = str(spec["reference_b"])
    n_a = int((subset["best_shape"] == ref_a).sum())
    n_b = int((subset["best_shape"] == ref_b).sum())
    ax.set_title(f"CN {cn}  (N = {len(subset):,})", loc="left", pad=5)
    ax.set_xlabel(rf"$\log_{{10}}[S(\mathrm{{{ref_a}}}) + 0.1]$")
    ax.set_ylabel(rf"$\log_{{10}}[S(\mathrm{{{ref_b}}}) + 0.1]$")
    ax.text(
        0.035,
        0.965,
        f"{ref_a}  {n_a:,}",
        transform=ax.transAxes,
        color=str(spec["color_a"]),
        fontsize=6.5,
        fontweight="bold",
        ha="left",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
        zorder=5,
    )
    ax.text(
        0.965,
        0.045,
        f"{ref_b}  {n_b:,}",
        transform=ax.transAxes,
        color=str(spec["color_b"]),
        fontsize=6.5,
        fontweight="bold",
        ha="right",
        va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
        zorder=5,
    )
    panel_label(ax, letter)
    return density


def draw_boundary_panel(ax: mpl.axes.Axes, summary: pd.DataFrame) -> None:
    panel_label(ax, "D")
    ax.set_title("Near-boundary records", loc="left", pad=5)
    order = ["CN4", "CN5", "CN6", "CN4-CN6"]
    labels = ["CN 4", "CN 5", "CN 6", "All"]
    colors = [NAVY, TEAL, ORANGE, INK]
    y = np.arange(len(order))[::-1]
    rows = summary.set_index("scope").loc[order]
    values = rows["boundary_margin_lt_1_percent"].to_numpy(dtype=float)
    numerators = rows["boundary_margin_lt_1_n"].to_numpy(dtype=int)
    denominators = rows["n"].to_numpy(dtype=int)
    ax.hlines(y, 0, values, color=GRID, linewidth=2.0, zorder=1)
    for yi, value, numerator, denominator, color in zip(
        y, values, numerators, denominators, colors
    ):
        ax.scatter(
            [value],
            [yi],
            s=28,
            facecolor=color,
            edgecolor="white",
            linewidth=0.55,
            zorder=3,
        )
        if value > 8:
            text_x = value - 0.30
            horizontal_alignment = "right"
        else:
            text_x = value + 0.35
            horizontal_alignment = "left"
        ax.text(
            text_x,
            yi,
            f"{value:.2f}%  ({numerator:,}/{denominator:,})",
            fontsize=6.6,
            ha=horizontal_alignment,
            va="center",
            color=INK,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.0},
        )
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 12.8)
    ax.set_ylim(-0.65, 3.65)
    ax.set_xticks([0, 3, 6, 9, 12])
    ax.set_xlabel(r"Records with $\Delta S < 1$ (%)")
    ax.text(
        0.02,
        0.04,
        r"$\Delta S=S_{\mathrm{second}}-S_{\mathrm{best}}$ (unrounded)",
        transform=ax.transAxes,
        fontsize=6.3,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    clean_axis(ax, grid_axis="x")


def write_figure(panel_data: pd.DataFrame, summary: pd.DataFrame) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(178 / 25.4, 129 / 25.4))
    fig.subplots_adjust(
        left=0.095,
        right=0.900,
        top=0.935,
        bottom=0.115,
        wspace=0.34,
        hspace=0.42,
    )
    density = draw_shape_panel(axes[0, 0], panel_data, 4, "A")
    draw_shape_panel(axes[0, 1], panel_data, 5, "B")
    draw_shape_panel(axes[1, 0], panel_data, 6, "C")
    draw_boundary_panel(axes[1, 1], summary)

    color_axis = fig.add_axes([0.919, 0.579, 0.014, 0.305])
    colorbar = fig.colorbar(density, cax=color_axis)
    colorbar.set_label("Records per hexbin", fontsize=6.7, labelpad=3)
    colorbar.ax.tick_params(labelsize=6.2, width=0.55, length=2.5)
    colorbar.outline.set_linewidth(0.55)
    fig.text(
        0.095,
        0.035,
        "Fixed tmQM/GFN2-xTB cohort: 48,057 records; 35,485 supported CN2-CN6 records; A-C show 33,863 CN4-CN6 records.",
        fontsize=6.25,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        ("_600dpi.png", {"dpi": 600}),
    ):
        fig.savefig(OUT / f"{STEM}{suffix}", facecolor="white", **kwargs)
    plt.close(fig)

    png = OUT / f"{STEM}_600dpi.png"
    with Image.open(png) as image:
        rgb = image.convert("RGB")
        rgb.save(png, dpi=(600, 600), compress_level=9)


def write_manifest() -> None:
    output_names = [
        "Figure3_caption.txt",
        f"{STEM}.pdf",
        f"{STEM}.svg",
        f"{STEM}_600dpi.png",
        "Figure3_panel_source.csv",
        "Figure3_summary_and_internal_QC.csv",
        "README.md",
    ]
    output_paths = [OUT / name for name in output_names]
    with (OUT / "Figure3_SHA256SUMS.txt").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        for path in output_paths:
            handle.write(f"{sha256(path)}  {path.name}\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with SOURCE_SUMMARY.open(encoding="utf-8") as handle:
        summary_json = json.load(handle)
    frame = pd.read_csv(SOURCE)
    validate_source(frame, summary_json)
    panel_data = write_panel_source(frame)
    summary = build_summary(panel_data)
    write_figure(panel_data, summary)
    write_manifest()
    print(f"Wrote corrected Figure 3 assets to {OUT}")


if __name__ == "__main__":
    main()
