#!/usr/bin/env python3
"""
Supplementary Figure S6 — Masked-field learning diagnostics and donor-field attribution.

Input files:
  fig4_field_learning_revision/fig4B_syntax_validity_vs_steps.csv
  fig4_field_learning_revision/fig4B_mask_ratio_recovery.csv
  fig4_field_learning_revision/fig4C_donor_field_attribution.csv
  masked_field_learning/donor_recovery_by_cn.csv
  masked_field_learning/donor_recovery_by_element.csv

Output: revision_results/si_figures/Fig_S6.{pdf,svg,png}

Panels:
  A — Training dynamics: syntax validity vs steps
  B — Mask-ratio recovery (shuffled-context baseline)
  C — Donor-field attribution controls (bar chart)
  D — Full-context donor-marker recovery by coordination number
  E — Full-context donor-marker recovery by donor element
"""
import sys, csv, textwrap
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "si_figures"))

from jacs_style import (apply_jacs_style, save_fig, label_panel, clean_axes,
                         light_grid, BLACK, DARK_GREY, MID_GREY, BLUE_GREY,
                         ACCENT, PALE_GREY, FAINT_GREY, SERIES_3, SERIES_4)
import matplotlib.pyplot as plt

apply_jacs_style()

FIG4   = ROOT / "revision_results" / "fig4_field_learning_revision"
# Canonical path references (actual files still in tool_a_ablation/)
MFL    = ROOT / "revision_results" / "tool_a_ablation"
OUTDIR = ROOT / "revision_results" / "si_figures"

def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

# ── Load data ────────────────────────────────────────
syntax = load_csv(FIG4 / "fig4B_syntax_validity_vs_steps.csv")
mask_ratio = load_csv(FIG4 / "fig4B_mask_ratio_recovery.csv")
donor_attr = load_csv(FIG4 / "fig4C_donor_field_attribution.csv")
by_cn = load_csv(MFL / "tool_a_ablation_by_cn.csv")
by_elem = load_csv(MFL / "tool_a_ablation_by_donor_element.csv")

# ── Figure layout: 3 rows ────────────────────────────
fig = plt.figure(figsize=(7.0, 8.0))
gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.35,
                      left=0.10, right=0.96, top=0.96, bottom=0.06)

# ═══════════════════════════════════════════════════════
# Panel A — Syntax validity vs training steps
# ═══════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0, 0])
label_panel(ax_a, "A")
clean_axes(ax_a)
light_grid(ax_a)

steps = [int(r["step"]) for r in syntax]
parse_v = [float(r["parse_validity_percent"]) for r in syntax]
strict_v = [float(r["strict_validity_percent"]) for r in syntax]

ax_a.plot(steps, parse_v, color=DARK_GREY, label="Parse valid", linewidth=1.3)
ax_a.plot(steps, strict_v, color=BLUE_GREY, label="Strict valid", linewidth=1.3)
ax_a.set_xlabel("Training step")
ax_a.set_ylabel("Validity (%)")
ax_a.set_title("Syntax learnability", fontsize=8)
ax_a.legend(loc="lower right", fontsize=6.5)
ax_a.set_ylim(-2, 105)

# ═══════════════════════════════════════════════════════
# Panel B — Mask-ratio recovery
# ═══════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[0, 1])
label_panel(ax_b, "B")
clean_axes(ax_b)
light_grid(ax_b)

ratios = [float(r["mask_ratio"]) for r in mask_ratio]
top1 = [float(r["top1_accuracy"]) for r in mask_ratio]
top5 = [float(r["top5_accuracy"]) for r in mask_ratio]
rand_bl = [float(r["random_baseline"]) for r in mask_ratio]

ax_b.plot(ratios, top1, marker="o", markersize=3.5, color=DARK_GREY,
          label="Top-1", linewidth=1.3)
ax_b.plot(ratios, top5, marker="s", markersize=3.5, color=BLUE_GREY,
          label="Top-5", linewidth=1.3)
ax_b.plot(ratios, rand_bl, ls="--", color=PALE_GREY,
          label="Shuffled context", linewidth=0.8)
ax_b.set_xlabel("Mask ratio")
ax_b.set_ylabel("Recovery (%)")
ax_b.set_title("Mask-ratio sensitivity", fontsize=8)
ax_b.legend(loc="upper right", fontsize=6)

# ═══════════════════════════════════════════════════════
# Panel C — Donor-field attribution controls
# ═══════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[1, 0])
label_panel(ax_c, "C")
clean_axes(ax_c)
light_grid(ax_c)

cond_order = [
    "full_context_MLM",
    "LigandFreq_lookup",
    "no_geometry_stereo",
    "ligand_SMILES_masked",
    "shuffled_ligand_SMILES",
]
cond_labels = [
    "Full context",
    "LigandFreq",
    "No geom/stereo",
    "Ligand masked",
    "Shuffled",
]
attr_map = {r["condition"]: r for r in donor_attr}
vals_c = [float(attr_map[c]["donor_top1"]) for c in cond_order]

x = np.arange(len(cond_labels))
# Grey/blue-grey for most; orange only for the key "Full context" bar
colors_c = [ACCENT] + [DARK_GREY, DARK_GREY, BLUE_GREY, BLUE_GREY]
bars = ax_c.bar(x, vals_c, width=0.55, color=colors_c, edgecolor=BLACK, linewidth=0.4)
ax_c.set_xticks(x)
ax_c.set_xticklabels(cond_labels, fontsize=5.5, rotation=25, ha="right")
ax_c.set_ylabel("Donor-marker Top-1 (%)")
ax_c.set_title("Donor-field attribution controls", fontsize=8)
for b, v in zip(bars, vals_c):
    ax_c.text(b.get_x() + b.get_width()/2, v + 0.8, f"{v:.1f}",
              ha="center", va="bottom", fontsize=5.5, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel D — Full-context donor-marker recovery by CN
# ═══════════════════════════════════════════════════════
ax_d = fig.add_subplot(gs[1, 1])
label_panel(ax_d, "D")
clean_axes(ax_d)
light_grid(ax_d)

cn_full = [r for r in by_cn if r["ablation_mode"] == "full_context"
           and int(r["CN"]) in (2, 3, 4, 5, 6) and int(r["n"]) >= 50]
cn_full.sort(key=lambda r: int(r["CN"]))

cn_labels_d = [f"CN={r['CN']}" for r in cn_full]
cn_vals = [float(r["top1_acc"]) * 100 for r in cn_full]
cn_n = [int(r["n"]) for r in cn_full]

x = np.arange(len(cn_labels_d))
bars_d = ax_d.bar(x, cn_vals, width=0.55, color=DARK_GREY, edgecolor=BLACK, linewidth=0.4)
ax_d.set_xticks(x)
ax_d.set_xticklabels(cn_labels_d, fontsize=6.5)
ax_d.set_ylabel("Donor-marker Top-1 (%)")
ax_d.set_title("Full-context donor recovery by CN", fontsize=8)
for b, v, n in zip(bars_d, cn_vals, cn_n):
    ax_d.text(b.get_x() + b.get_width()/2, v + 0.5, f"{v:.1f}\n(n={n:,})",
              ha="center", va="bottom", fontsize=5, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel E — Full-context donor-marker recovery by element
# ═══════════════════════════════════════════════════════
ax_e = fig.add_subplot(gs[2, :])
label_panel(ax_e, "E", x=-0.05)
clean_axes(ax_e)
light_grid(ax_e)

elem_full = [r for r in by_elem if r["ablation_mode"] == "full_context"
             and int(r["n"]) >= 20]
elem_full.sort(key=lambda r: -float(r["top1_acc"]))

elem_labels = [r["donor_element"] for r in elem_full]
elem_vals = [float(r["top1_acc"]) * 100 for r in elem_full]
elem_n = [int(r["n"]) for r in elem_full]

# Weighted average for verification
wavg_e = sum(v * n for v, n in zip(
    [float(r["top1_acc"]) for r in elem_full], elem_n)) / sum(elem_n) * 100

x = np.arange(len(elem_labels))
# Grey/blue-grey for all bars; no orange highlight
bars_e = ax_e.bar(x, elem_vals, width=0.6, color=DARK_GREY, edgecolor=BLACK, linewidth=0.4)
ax_e.set_xticks(x)
ax_e.set_xticklabels(elem_labels, fontsize=6.5)
ax_e.set_ylabel("Donor-marker Top-1 (%)")
ax_e.set_title("Full-context donor-marker recovery by donor element", fontsize=8)
ax_e.set_ylim(0, 105)
for b, v, n in zip(bars_e, elem_vals, elem_n):
    ax_e.text(b.get_x() + b.get_width()/2, v + 0.5, f"n={n:,}",
              ha="center", va="bottom", fontsize=4.5, color=MID_GREY, rotation=45)
# Annotate weighted average
ax_e.axhline(wavg_e, ls=":", lw=0.7, color=ACCENT, zorder=0)
ax_e.text(len(elem_labels) - 0.5, wavg_e + 1.5,
          f"weighted avg = {wavg_e:.1f}%",
          ha="right", fontsize=6, color=ACCENT)

# ── Save ─────────────────────────────────────────────
save_fig(fig, "Fig_S6", OUTDIR)
plt.close(fig)

# ── Notes ────────────────────────────────────────────
notes = textwrap.dedent(f"""\
    Fig_S6 — Masked-field learning diagnostics and donor-field attribution

    Panel A: Syntax validity vs training steps
      Source: fig4_field_learning_revision/fig4B_syntax_validity_vs_steps.csv
      x = training step; y = validity (%)
      Two lines: parse valid (dark grey), strict valid (blue-grey)

    Panel B: Mask-ratio recovery
      Source: fig4_field_learning_revision/fig4B_mask_ratio_recovery.csv
      x = mask ratio; y = recovery (%)
      Lines: Top-1 (dark grey), Top-5 (blue-grey)
      Dashed line: shuffled-context baseline (21.6%) — ligand SMILES shuffled
      to random positions, confirming near-chance performance without correct
      ligand context. This is NOT a uniform-random-over-vocabulary baseline.

    Panel C: Donor-field attribution controls
      Source: fig4_field_learning_revision/fig4C_donor_field_attribution.csv
      Bar chart: 5 conditions; y = donor-marker Top-1 accuracy (%)
      Full context = {vals_c[0]:.1f}% (n = 19,992); orange highlight

    Panel D: Full-context donor-marker recovery by coordination number
      Source: masked_field_learning/donor_recovery_by_cn.csv
      (actual file: tool_a_ablation/tool_a_ablation_by_cn.csv, full_context rows)
      CN ∈ {{2,3,4,5,6}} with n ≥ 50
      Metric: donor-marker Top-1, same evaluation protocol as Panel C

    Panel E: Full-context donor-marker recovery by donor element
      Source: masked_field_learning/donor_recovery_by_element.csv
      (actual file: tool_a_ablation/tool_a_ablation_by_donor_element.csv, full_context rows)
      Same metric and evaluation set as Panel C (n = 19,992, weighted avg = {wavg_e:.1f}%)
      Bars sorted by descending accuracy; n ≥ 20 filter
      Orange horizontal line = weighted average across all elements
""")
(OUTDIR / "Fig_S6_notes.txt").write_text(notes)
print("Done: Fig_S6")
