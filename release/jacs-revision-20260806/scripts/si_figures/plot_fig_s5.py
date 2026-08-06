#!/usr/bin/env python3
"""
Supplementary Figure S5 — Full-CSD boundary records and L3-family statistics.

Input files:
  csd_pathfinder_full/full_csd_pathfinder_summary.json
  csd_pathfinder_full/l3_family_geometry_trajectories.csv
  csd_pathfinder_full/cn4_sp_td_atlas.csv
  csd_pathfinder_full/cn5_tbpy_spy_atlas.csv
  csd_pathfinder_full/cn6_oh_distortion_atlas.csv

Output: revision_results/si_figures/Fig_S5.{pdf,svg,png}

Panels:
  A — Boundary fraction by CN (bar, CN 4/5/6 only)
  B — Dominant boundary shape pairs (horizontal bar, grey/blue-grey)
  C — L3 family size distribution (histogram, log x)
  D — L1-state diversity within L3 families (histogram)
  E — [CuCl₄]²⁻ L3 family: 1-D geometry-coordinate strip + zoom inset
"""
import sys, csv, json, textwrap
from pathlib import Path
from collections import Counter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "si_figures"))

from jacs_style import (apply_jacs_style, save_fig, label_panel, clean_axes,
                         light_grid, BLACK, DARK_GREY, MID_GREY, BLUE_GREY,
                         ACCENT, PALE_GREY, FAINT_GREY)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

apply_jacs_style()

CSD    = ROOT / "revision_results" / "csd_pathfinder_full"
OUTDIR = ROOT / "revision_results" / "si_figures"

# ── Load summary ─────────────────────────────────────
with open(CSD / "full_csd_pathfinder_summary.json") as f:
    summary = json.load(f)

cn_dist = summary["cn_distribution"]  # str keys

# ── Load atlases for boundary computation ────────────
def load_atlas_full(path, col_x, col_y):
    xs, ys, refs, l3s = [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                xs.append(float(row[col_x]))
                ys.append(float(row[col_y]))
                refs.append(row["refcode"])
                l3s.append(row.get("L3_hash", ""))
            except (ValueError, KeyError):
                pass
    return np.array(xs), np.array(ys), refs, l3s

cn4_x, cn4_y, cn4_refs, cn4_l3 = load_atlas_full(
    CSD / "cn4_sp_td_atlas.csv", "S_Td", "S_SP")
cn5_x, cn5_y, _, _ = load_atlas_full(
    CSD / "cn5_tbpy_spy_atlas.csv", "S_TBPY", "S_SPY")
cn6_x, cn6_y, _, _ = load_atlas_full(
    CSD / "cn6_oh_distortion_atlas.csv", "S_TPr", "S_Oh")

d4 = np.abs(cn4_x - cn4_y)
d5 = np.abs(cn5_x - cn5_y)
d6 = np.abs(cn6_x - cn6_y)

# ── Load L3 families ─────────────────────────────────
fam = []
with open(CSD / "l3_family_geometry_trajectories.csv") as f:
    for row in csv.DictReader(f):
        fam.append(row)

# ── Figure layout: 3 rows × 2 cols ──────────────────
fig = plt.figure(figsize=(7.0, 8.5))
gs = fig.add_gridspec(3, 2, hspace=0.48, wspace=0.35,
                      left=0.10, right=0.96, top=0.96, bottom=0.05)

# ═══════════════════════════════════════════════════════
# Panel A — Boundary fraction by CN (4/5/6 only)
# ═══════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0, 0])
label_panel(ax_a, "A")
clean_axes(ax_a)
light_grid(ax_a)

cn_labels = ["4", "5", "6"]
bf = {
    "4": np.sum(d4 < 1.0) / len(d4) * 100 if len(d4) else 0,
    "5": np.sum(d5 < 1.0) / len(d5) * 100 if len(d5) else 0,
    "6": np.sum(d6 < 1.0) / len(d6) * 100 if len(d6) else 0,
}
n_vals = {"4": len(d4), "5": len(d5), "6": len(d6)}
fracs = [bf[c] for c in cn_labels]
bar_colors_a = [DARK_GREY, MID_GREY, BLUE_GREY]
x = np.arange(len(cn_labels))
bars = ax_a.bar(x, fracs, width=0.55, color=bar_colors_a,
                edgecolor=BLACK, linewidth=0.4)
ax_a.set_xticks(x)
ax_a.set_xticklabels([f"CN = {c}\n(n={n_vals[c]:,})" for c in cn_labels], fontsize=6.5)
ax_a.set_ylabel("Boundary fraction (%)")
ax_a.set_title("Boundary records by CN", fontsize=8)
for b, v in zip(bars, fracs):
    ax_a.text(b.get_x() + b.get_width()/2, v + 0.5, f"{v:.1f}%",
              ha="center", va="bottom", fontsize=6, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel B — Dominant boundary shape pairs
# ═══════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[0, 1])
label_panel(ax_b, "B")
clean_axes(ax_b)
light_grid(ax_b, axis="x")

pairs = [
    ("SP \u2194 Td",       int(np.sum(d4 < 1.0))),
    ("SPY \u2194 TBPY",    int(np.sum(d5 < 1.0))),
    ("Oh \u2194 TPr",      int(np.sum(d6 < 1.0))),
]
pair_labels = [p[0] for p in pairs]
pair_counts = [p[1] for p in pairs]
pair_colors = [DARK_GREY, MID_GREY, BLUE_GREY]
y_b = np.arange(len(pairs))
bars_b = ax_b.barh(y_b, pair_counts, height=0.50, color=pair_colors,
                    edgecolor=BLACK, linewidth=0.4)
ax_b.set_yticks(y_b)
ax_b.set_yticklabels(pair_labels, fontsize=7)
ax_b.set_xlabel("Boundary records (\u0394CShM < 1)")
ax_b.set_title("Dominant boundary shape pairs", fontsize=8)
for b, v in zip(bars_b, pair_counts):
    ax_b.text(v + max(pair_counts) * 0.02, b.get_y() + b.get_height()/2,
              f"{v:,}", va="center", fontsize=6, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel C — L3 family size distribution
# ═══════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[1, 0])
label_panel(ax_c, "C")
clean_axes(ax_c)
light_grid(ax_c)

sizes = [int(r["n_records"]) for r in fam if int(r["n_records"]) >= 2]
bins = np.logspace(np.log10(1.5), np.log10(max(sizes) + 1), 30)
ax_c.hist(sizes, bins=bins, color=PALE_GREY, edgecolor=DARK_GREY, linewidth=0.5)
ax_c.set_xscale("log")
ax_c.set_xlabel("Records per L3 family")
ax_c.set_ylabel("Number of families")
ax_c.set_title("L3 family size distribution", fontsize=8)
ax_c.text(0.97, 0.92, f"{len(sizes):,} nontrivial\nL3 families",
          transform=ax_c.transAxes, ha="right", va="top",
          fontsize=6, color=MID_GREY)

# ═══════════════════════════════════════════════════════
# Panel D — L1-state diversity within L3 families
# ═══════════════════════════════════════════════════════
ax_d = fig.add_subplot(gs[1, 1])
label_panel(ax_d, "D")
clean_axes(ax_d)
light_grid(ax_d)

n_l1 = [int(r["unique_L1"]) for r in fam if int(r["n_records"]) >= 2]
max_l1 = max(n_l1) if n_l1 else 1
bins_d = np.arange(0.5, min(max_l1 + 1.5, 15.5), 1)
ax_d.hist(n_l1, bins=bins_d, color=BLUE_GREY, edgecolor=BLACK, linewidth=0.4)
ax_d.set_xlabel("Unique L1 states per L3 family")
ax_d.set_ylabel("Number of families")
ax_d.set_title("L1-state diversity within L3 families", fontsize=8)
multi_l1 = sum(1 for v in n_l1 if v >= 2)
ax_d.text(0.97, 0.92,
          f"{multi_l1:,} families with\n\u2265 2 geometry states",
          transform=ax_d.transAxes, ha="right", va="top",
          fontsize=6, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel E — [CuCl₄]²⁻ L3 family: scatter + zoom inset
# ═══════════════════════════════════════════════════════
ax_e = fig.add_subplot(gs[2, :])
label_panel(ax_e, "E", x=-0.05)
clean_axes(ax_e)

# Find CuCl4 family (largest CN=4 boundary-crossing family)
cn4_fam = [r for r in fam if r["CN"] == "4" and r["crosses_shape_boundary"] == "True"
           and int(r["n_records"]) >= 3]
cn4_fam.sort(key=lambda r: -int(r["n_records"]))
chosen = cn4_fam[0]
l3_hash = chosen["L3_hash"]

# Collect all family members from atlas
fam_mask = np.array([h == l3_hash for h in cn4_l3])
fam_td = cn4_x[fam_mask]
fam_sp = cn4_y[fam_mask]
fam_delta = np.abs(fam_td - fam_sp)
fam_refs_arr = [cn4_refs[i] for i in range(len(cn4_refs)) if fam_mask[i]]
n_boundary = int(np.sum(fam_delta < 1.0))
n_total = len(fam_td)

# Background: faint scatter of all CN4 (clipped to 99th pct)
bg_mask = (cn4_x < 22) & (cn4_y < 20)
ax_e.scatter(cn4_x[bg_mask], cn4_y[bg_mask], s=0.8, c=FAINT_GREY, alpha=0.25,
             zorder=1, rasterized=True)
ax_e.plot([0, 22], [0, 22], ls="--", lw=0.5, color=PALE_GREY, zorder=2)

# Family members — colour by boundary status
is_boundary = fam_delta < 1.0
ax_e.scatter(fam_td[~is_boundary], fam_sp[~is_boundary], s=12, c=MID_GREY,
             edgecolors=DARK_GREY, linewidths=0.3, zorder=3, label="Clear")
ax_e.scatter(fam_td[is_boundary], fam_sp[is_boundary], s=18, c=ACCENT,
             edgecolors=BLACK, linewidths=0.4, zorder=4, label="Boundary")

ax_e.set_xlim(0, 5)
ax_e.set_ylim(0, 12)
ax_e.set_xlabel("CShM(Td)")
ax_e.set_ylabel("CShM(SP)")
ax_e.set_title(
    f"[CuCl$_4$]$^{{2\u2212}}$ L3 family: {n_total} members, "
    f"{n_boundary} at boundary",
    fontsize=7.5, pad=3)
ax_e.legend(loc="upper left", fontsize=6, markerscale=1.5,
            framealpha=0.8, edgecolor=PALE_GREY)

# Zoom inset on the boundary region
ax_ins = inset_axes(ax_e, width="38%", height="50%", loc="center right",
                    borderpad=1.5)
clean_axes(ax_ins)
# Boundary zone: ΔCShM < 2 region
zoom_mask = fam_delta < 2.0
ax_ins.scatter(fam_td[zoom_mask], fam_sp[zoom_mask], s=25,
               c=[ACCENT if d < 1.0 else MID_GREY for d in fam_delta[zoom_mask]],
               edgecolors=BLACK, linewidths=0.3, zorder=4)
# Label a few refcodes in zoom
zoom_refs = [fam_refs_arr[i] for i in range(len(fam_refs_arr))
             if fam_delta[i] < 2.0]
zoom_td = fam_td[zoom_mask]
zoom_sp = fam_sp[zoom_mask]
# Label every 8th point to avoid clutter
for idx in range(0, len(zoom_refs), max(1, len(zoom_refs)//6)):
    ax_ins.annotate(zoom_refs[idx], (zoom_td[idx], zoom_sp[idx]),
                    fontsize=4, color=DARK_GREY,
                    xytext=(3, 3), textcoords="offset points")
# Diagonal in inset
zmax = max(zoom_td.max(), zoom_sp.max()) + 0.5
ax_ins.plot([0, zmax], [0, zmax], ls="--", lw=0.4, color=PALE_GREY, zorder=2)
ax_ins.set_xlabel("CShM(Td)", fontsize=5.5)
ax_ins.set_ylabel("CShM(SP)", fontsize=5.5)
ax_ins.tick_params(labelsize=5)
ax_ins.set_title("Boundary zoom", fontsize=6, pad=2)
# Indicate zoom region on main
from matplotlib.patches import FancyBboxPatch
rect = plt.Rectangle((zoom_td.min() - 0.1, zoom_sp.min() - 0.2),
                       zoom_td.max() - zoom_td.min() + 0.2,
                       zoom_sp.max() - zoom_sp.min() + 0.4,
                       fill=False, edgecolor=ACCENT, lw=0.7, ls="--", zorder=5)
ax_e.add_patch(rect)

# ── Save ─────────────────────────────────────────────
save_fig(fig, "Fig_S5", OUTDIR)
plt.close(fig)

# ── Notes ────────────────────────────────────────────
notes = textwrap.dedent(f"""\
    Fig_S5 — Full-CSD boundary records and L3-family statistics

    Panel A: Boundary fraction by CN
      Source: cn4/5/6 atlas CSVs (\u0394CShM computed); threshold = 1.0
      Only CN = 4, 5, 6 shown (CN 2/3 not evaluated: no dual-shape competition)

    Panel B: Dominant boundary shape pairs
      Source: same atlas CSVs; boundary count per shape pair
      Horizontal bar chart; count of entries with \u0394CShM < 1.0
      Only the three dominant channels shown (one per evaluated CN)

    Panel C: L3 family size distribution
      Source: l3_family_geometry_trajectories.csv
      Families with \u2265 2 members; log-x histogram
      {len(sizes):,} nontrivial L3 families

    Panel D: L1-state diversity within L3 families
      Source: same CSV; unique_L1 per family
      Integer histogram

    Panel E: [CuCl\u2084]\u00b2\u207b L3 family — geometry spread
      Source: cn4_sp_td_atlas.csv; L3_hash = {l3_hash}
      {n_total} members, {n_boundary} at boundary (\u0394CShM < 1.0)
      Main scatter: all CN = 4 background (faint) + family members
      Zoom inset: boundary region with refcode labels
      Note: Geometry spread is a static comparison among independently observed
      crystal structures; it does not imply a dynamic interconversion pathway.
""")
(OUTDIR / "Fig_S5_notes.txt").write_text(notes)
print("Done: Fig_S5")
