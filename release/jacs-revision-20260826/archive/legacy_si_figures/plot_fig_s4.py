#!/usr/bin/env python3
"""
Supplementary Figure S4 — Extended CShM geometry atlases and threshold sensitivity.

Input files:
  csd_pathfinder_full/cn4_sp_td_atlas.csv
  csd_pathfinder_full/cn5_tbpy_spy_atlas.csv
  csd_pathfinder_full/cn6_oh_distortion_atlas.csv
  csd_pathfinder_full/full_csd_pathfinder_summary.json

Output: revision_results/si_figures/Fig_S4.{pdf,svg,png}

Panels:
  A — CN=4 SP–Td hexbin atlas
  B — CN=5 SPY–TBPY hexbin atlas
  C — CN=6 Oh–TPr hexbin atlas
  D — Boundary-threshold sensitivity (line plot)
"""
import sys, csv, json, textwrap
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "si_figures"))

from jacs_style import (apply_jacs_style, save_fig, label_panel, clean_axes,
                         light_grid, BLACK, DARK_GREY, MID_GREY, BLUE_GREY,
                         ACCENT, PALE_GREY, FAINT_GREY, ACCENT_LT)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

apply_jacs_style()

CSD    = ROOT / "revision_results" / "csd_pathfinder_full"
OUTDIR = ROOT / "revision_results" / "si_figures"

# ── Load atlases ─────────────────────────────────────
def load_atlas(path, col_x, col_y):
    x, y = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                x.append(float(row[col_x]))
                y.append(float(row[col_y]))
            except (ValueError, KeyError):
                continue
    return np.array(x), np.array(y)

cn4_x, cn4_y = load_atlas(CSD / "cn4_sp_td_atlas.csv", "S_Td", "S_SP")
cn5_x, cn5_y = load_atlas(CSD / "cn5_tbpy_spy_atlas.csv", "S_TBPY", "S_SPY")
cn6_x, cn6_y = load_atlas(CSD / "cn6_oh_distortion_atlas.csv", "S_TPr", "S_Oh")

print(f"  CN4: {len(cn4_x)} pts, CN5: {len(cn5_x)} pts, CN6: {len(cn6_x)} pts")

# ── Restrained grey-to-blue-grey colormap (no rainbow) ─
cmap_colors = ["#FFFFFF", "#E5E7EA", "#B8BFC8", "#8A95A3",
               "#5C7088", "#3E5266"]
cmap = mcolors.LinearSegmentedColormap.from_list("jacs_grey_blue", cmap_colors, N=256)

# ── Axis limits (approximately 99th percentile, rounded readably) ──
AXIS_LIMITS = {
    "CN4": (0, 22, 0, 20),   # x_lo, x_hi, y_lo, y_hi
    "CN5": (0, 13, 0, 16),
    "CN6": (3, 14, 0, 8),
}

# ── Figure layout ────────────────────────────────────
fig = plt.figure(figsize=(7.0, 7.0))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34,
                      left=0.09, right=0.96, top=0.96, bottom=0.07)

# ═══════════════════════════════════════════════════════
# Helper: hexbin atlas panel with log(count+1) density
# ═══════════════════════════════════════════════════════
def plot_atlas(ax, x, y, xlabel, ylabel, label_lo, label_hi,
               panel_label, cn_label, xlim, ylim):
    label_panel(ax, panel_label)
    clean_axes(ax)
    # Clip to display range for hex binning
    mask = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    xc, yc = x[mask], y[mask]
    hb = ax.hexbin(xc, yc, gridsize=50, cmap=cmap, mincnt=1,
                   linewidths=0.1, edgecolors="none",
                   bins="log")  # log(count+1) colour scaling
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(f"CN = {cn_label}  (n = {len(xc):,})", fontsize=8, pad=4)
    # Thin light-grey diagonal
    diag_max = max(xlim[1], ylim[1])
    ax.plot([0, diag_max], [0, diag_max], ls="--", lw=0.5, color=PALE_GREY, zorder=5)
    # Region labels
    ax.text(0.05, 0.93, label_lo, transform=ax.transAxes,
            fontsize=6, color=DARK_GREY, fontstyle="italic", va="top")
    ax.text(0.93, 0.05, label_hi, transform=ax.transAxes,
            fontsize=6, color=DARK_GREY, fontstyle="italic", ha="right")
    return hb

ax_a = fig.add_subplot(gs[0, 0])
hb_a = plot_atlas(ax_a, cn4_x, cn4_y,
           "CShM(Td)", "CShM(SP)", "Td-like", "SP-like",
           "A", "4", AXIS_LIMITS["CN4"][:2], AXIS_LIMITS["CN4"][2:])

ax_b = fig.add_subplot(gs[0, 1])
hb_b = plot_atlas(ax_b, cn5_x, cn5_y,
           "CShM(TBPY)", "CShM(SPY)", "TBPY-like", "SPY-like",
           "B", "5", AXIS_LIMITS["CN5"][:2], AXIS_LIMITS["CN5"][2:])

ax_c = fig.add_subplot(gs[1, 0])
hb_c = plot_atlas(ax_c, cn6_x, cn6_y,
           "CShM(TPr)", "CShM(Oh)", "TPr-like", "Oh-like",
           "C", "6", AXIS_LIMITS["CN6"][:2], AXIS_LIMITS["CN6"][2:])

# Shared colorbar — labelled as log(count+1) relative density
cbar_ax = fig.add_axes([0.10, 0.025, 0.36, 0.012])
cb = fig.colorbar(hb_c, cax=cbar_ax, orientation="horizontal")
cb.set_label("log$_{10}$(count + 1)", fontsize=6)
cb.ax.tick_params(labelsize=5.5)

# ═══════════════════════════════════════════════════════
# Panel D — Boundary-threshold sensitivity
# ═══════════════════════════════════════════════════════
ax_d = fig.add_subplot(gs[1, 1])
label_panel(ax_d, "D")
clean_axes(ax_d)
light_grid(ax_d)

# Compute boundary fractions at different thresholds from atlas data
thresholds = np.arange(0.0, 3.1, 0.25)

def boundary_frac(x, y, thr):
    delta = np.abs(x - y)
    return np.sum(delta < thr) / len(delta) * 100 if len(delta) > 0 else 0

# Muted grey / blue-grey lines for CN curves
cn_line_colors = [DARK_GREY, MID_GREY, BLUE_GREY]
for (data, lbl), color in zip([
    ((cn4_x, cn4_y), "CN = 4"),
    ((cn5_x, cn5_y), "CN = 5"),
    ((cn6_x, cn6_y), "CN = 6"),
], cn_line_colors):
    fracs = [boundary_frac(data[0], data[1], t) for t in thresholds]
    ax_d.plot(thresholds, fracs, marker="o", markersize=3, color=color,
              label=lbl, linewidth=1.3)

# Default threshold — muted orange accent
ax_d.axvline(1.0, ls=":", lw=1.0, color=ACCENT, zorder=0)
ax_d.text(1.08, ax_d.get_ylim()[1] * 0.92, "Δ = 1.0\n(default)", fontsize=5.5,
          color=ACCENT, va="top")
ax_d.set_xlabel("ΔCShM threshold")
ax_d.set_ylabel("Boundary fraction (%)")
ax_d.set_title("Threshold sensitivity", fontsize=8)
ax_d.legend(loc="upper left", fontsize=6.5)

# ── Save ─────────────────────────────────────────────
save_fig(fig, "Fig_S4", OUTDIR)
plt.close(fig)

# ── Notes ────────────────────────────────────────────
notes = textwrap.dedent(f"""\
    Fig_S4 — Extended CShM geometry atlases and threshold sensitivity

    Panel A: CN=4 SP–Td atlas
      Source: csd_pathfinder_full/cn4_sp_td_atlas.csv
      x = CShM(Td), y = CShM(SP); hexbin density, log(count+1) colour scaling
      Axis range clipped to ≈99th percentile; n = {len(cn4_x):,} entries
      Dashed diagonal = equal CShM (boundary-like region)

    Panel B: CN=5 SPY–TBPY atlas
      Source: csd_pathfinder_full/cn5_tbpy_spy_atlas.csv
      x = CShM(TBPY), y = CShM(SPY); hexbin density, log(count+1)
      n = {len(cn5_x):,} entries

    Panel C: CN=6 Oh–TPr atlas
      Source: csd_pathfinder_full/cn6_oh_distortion_atlas.csv
      x = CShM(TPr), y = CShM(Oh); hexbin density, log(count+1)
      n = {len(cn6_x):,} entries

    Shared colorbar: log₁₀(count+1) relative density; grey-to-blue-grey cmap.
    Each panel uses log-scaled hexbin counts; colorbar placed below Panel C.

    Panel D: Boundary-threshold sensitivity
      Source: computed from CN4/5/6 atlas data
      x = ΔCShM threshold; y = boundary fraction (%)
      Lines: CN 4/5/6 in muted grey/blue-grey tones
      Vertical dotted orange line at Δ=1.0 (default CoordRep threshold)
""")
(OUTDIR / "Fig_S4_notes.txt").write_text(notes)
print("Done: Fig_S4")
