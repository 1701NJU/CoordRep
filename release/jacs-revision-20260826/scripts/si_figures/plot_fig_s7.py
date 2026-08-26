#!/usr/bin/env python3
"""
Supplementary Figure S7 — Multidentate donor-set evaluation.

Input files:
  multidentate/donor_set_recovery_by_denticity.csv
  (actual file: multidentate/tool_a_by_denticity_final.csv)

Output: revision_results/si_figures/Fig_S7.{pdf,svg,png}

Panels:
  A — Metric schematic (token-level vs ligand-level exact match)
  B — Mono- and bidentate donor-set recovery (bar chart)
  C — Bidentate masking-regime comparison with shuffled-ligand control
"""
import sys, csv, json, textwrap
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "si_figures"))

from jacs_style import (apply_jacs_style, save_fig, label_panel, clean_axes,
                         light_grid, BLACK, DARK_GREY, MID_GREY, BLUE_GREY,
                         ACCENT, PALE_GREY, FAINT_GREY)
import matplotlib.pyplot as plt

apply_jacs_style()

MULTI  = ROOT / "revision_results" / "multidentate"
OUTDIR = ROOT / "revision_results" / "si_figures"

# ── Load data ────────────────────────────────────────
def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

dent = load_csv(MULTI / "tool_a_by_denticity_final.csv")

# ── Figure layout: 2 × 2 with Panel A spanning left, B right-top, C full bottom
fig = plt.figure(figsize=(7.0, 6.5))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.35,
                      left=0.09, right=0.97, top=0.96, bottom=0.08)

# ═══════════════════════════════════════════════════════
# Panel A — Metric schematic
# ═══════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_xlim(0, 100)
ax_a.set_ylim(0, 50)
ax_a.axis("off")
label_panel(ax_a, "A", x=-0.02, y=1.08)

ax_a.text(50, 48, "Token-level vs ligand-level evaluation",
          ha="center", va="top", fontsize=7.5, fontweight="bold", color=BLACK)

# Bidentate example
y0 = 39
ax_a.text(5, y0, "Bidentate ligand L1:", fontsize=6.5, color=DARK_GREY,
          fontweight="bold", va="center")
ax_a.text(5, y0 - 5, "True donor set:", fontsize=6.5, color=DARK_GREY, va="center")
ax_a.text(40, y0 - 5, "{L1:N:1, L1:N:2}", fontfamily="monospace",
          fontsize=7, color="#3A7A3A", va="center")

# Case 1: partial match → donor-set mismatch
ax_a.plot([5, 95], [y0 - 9, y0 - 9], color=PALE_GREY, lw=0.5)
ax_a.text(5, y0 - 12, "Predicted:", fontsize=6.5, color=DARK_GREY, va="center")
ax_a.text(40, y0 - 12, "{L1:N:1, L1:O:1}", fontfamily="monospace",
          fontsize=7, color=ACCENT, va="center")
ax_a.text(5, y0 - 16, "\u2192 Token-level:", fontsize=6, color=MID_GREY, va="center")
ax_a.text(40, y0 - 16, "1/2 correct (50%)", fontsize=6, color=MID_GREY, va="center")
ax_a.text(5, y0 - 19.5, "\u2192 Ligand-level:", fontsize=6, color=MID_GREY, va="center")
ax_a.text(40, y0 - 19.5, "no exact match (mismatch)", fontsize=6, color=ACCENT,
          fontweight="bold", va="center")

# Case 2: full match
ax_a.plot([5, 95], [y0 - 23, y0 - 23], color=PALE_GREY, lw=0.5)
ax_a.text(5, y0 - 26, "Predicted:", fontsize=6.5, color=DARK_GREY, va="center")
ax_a.text(40, y0 - 26, "{L1:N:1, L1:N:2}", fontfamily="monospace",
          fontsize=7, color="#3A7A3A", va="center")
ax_a.text(5, y0 - 30, "\u2192 Token-level:", fontsize=6, color=MID_GREY, va="center")
ax_a.text(40, y0 - 30, "2/2 correct (100%)", fontsize=6, color=MID_GREY, va="center")
ax_a.text(5, y0 - 33.5, "\u2192 Ligand-level:", fontsize=6, color=MID_GREY, va="center")
ax_a.text(40, y0 - 33.5, "exact match", fontsize=6, color="#3A7A3A",
          fontweight="bold", va="center")

# ═══════════════════════════════════════════════════════
# Panel B — Mono- and bidentate donor-set recovery
# ═══════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[0, 1])
label_panel(ax_b, "B")
clean_axes(ax_b)
light_grid(ax_b)

# full_context rows only (d=1 and d=2)
fc = [r for r in dent if r["ablation_mode"] == "full_context"]
fc.sort(key=lambda r: r["denticity"])
dent_labels = [r["denticity"].replace("dent=", "d = ") for r in fc]
exact_match = [float(r["per_ligand_all_correct"]) * 100 for r in fc]
per_donor = [float(r["per_donor_top1"]) * 100 for r in fc]

x = np.arange(len(dent_labels))
w = 0.35
bars_per = ax_b.bar(x - w/2, per_donor, width=w, color=PALE_GREY,
                     edgecolor=BLACK, linewidth=0.4, label="Per-donor Top-1")
bars_exact = ax_b.bar(x + w/2, exact_match, width=w, color=DARK_GREY,
                       edgecolor=BLACK, linewidth=0.4, label="Ligand exact match")
ax_b.set_xticks(x)
ax_b.set_xticklabels(dent_labels, fontsize=7)
ax_b.set_ylabel("Accuracy (%)")
ax_b.set_title("Mono- and bidentate donor-set recovery", fontsize=8)
ax_b.legend(loc="upper right", fontsize=5.5)
ax_b.set_ylim(0, 105)

for b, v in zip(bars_exact, exact_match):
    ax_b.text(b.get_x() + b.get_width()/2, v + 1, f"{v:.1f}",
              ha="center", va="bottom", fontsize=5.5, color=DARK_GREY)
for b, v in zip(bars_per, per_donor):
    ax_b.text(b.get_x() + b.get_width()/2, v + 1, f"{v:.1f}",
              ha="center", va="bottom", fontsize=5.5, color=MID_GREY)

# ═══════════════════════════════════════════════════════
# Panel C — Masking regime comparison (bidentate) with shuffled control
# ═══════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[1, :])
label_panel(ax_c, "C", x=-0.04)
clean_axes(ax_c)
light_grid(ax_c)

# All regimes for bidentate
regimes = {
    "full_context": "Full\ncontext",
    "no_geometry_stereo": "No geom/\nstereo",
    "no_ligand_smiles_keep_length": "Ligand\nmasked",
    "no_ligand_smiles_collapsed": "Ligand\ncollapsed",
    "metal_cn_only": "Metal+CN\nonly",
    "shuffled_ligand_smiles_control": "Shuffled\nligand",
}

bident = {}
for r in dent:
    if r["denticity"] == "dent=2" and r["ablation_mode"] in regimes:
        bident[r["ablation_mode"]] = {
            "per_donor": float(r["per_donor_top1"]) * 100,
            "exact": float(r["per_ligand_all_correct"]) * 100,
        }

regime_keys = [k for k in regimes if k in bident]
regime_labels = [regimes[k] for k in regime_keys]
per_donor_vals = [bident[k]["per_donor"] for k in regime_keys]
exact_vals = [bident[k]["exact"] for k in regime_keys]

x = np.arange(len(regime_labels))
w = 0.35
# Grey/blue-grey bars; orange only for full_context exact-match bar
pd_colors = [PALE_GREY] * len(regime_keys)
ex_colors = [ACCENT if k == "full_context" else BLUE_GREY for k in regime_keys]

ax_c.bar(x - w/2, per_donor_vals, width=w, color=pd_colors,
         edgecolor=BLACK, linewidth=0.4, label="Per-donor Top-1")
ax_c.bar(x + w/2, exact_vals, width=w, color=ex_colors,
         edgecolor=BLACK, linewidth=0.4, label="Ligand exact match")
ax_c.set_xticks(x)
ax_c.set_xticklabels(regime_labels, fontsize=6)
ax_c.set_ylabel("Accuracy (%)")
ax_c.set_title("Bidentate donor-set recovery by masking regime", fontsize=8)
ax_c.legend(loc="upper right", fontsize=6)
ax_c.set_ylim(0, 105)

for xi, (pd, ex) in enumerate(zip(per_donor_vals, exact_vals)):
    ax_c.text(xi - w/2, pd + 1, f"{pd:.1f}", ha="center", va="bottom",
              fontsize=5, color=DARK_GREY)
    ax_c.text(xi + w/2, ex + 1, f"{ex:.1f}", ha="center", va="bottom",
              fontsize=5, color=DARK_GREY)

# ── Save ─────────────────────────────────────────────
save_fig(fig, "Fig_S7", OUTDIR)
plt.close(fig)

# ── Notes ────────────────────────────────────────────
notes = textwrap.dedent("""\
    Fig_S7 — Multidentate donor-set evaluation

    Panel A: Metric schematic
      No source data — annotation-only panel
      Illustrates difference between token-level (per-donor) accuracy
      and ligand-level exact match for bidentate ligands.
      Partial recovery is labelled "no exact match (mismatch),"
      not "FAIL" or "invalid completion."

    Panel B: Mono- and bidentate donor-set recovery
      Source: multidentate/donor_set_recovery_by_denticity.csv
      (actual file: multidentate/tool_a_by_denticity_final.csv, full_context rows)
      Grouped bar: per-donor Top-1 (pale grey) vs ligand exact match (dark grey)
      x = denticity (d=1, d=2); y = accuracy (%)
      No tridentate or higher-denticity data available.

    Panel C: Bidentate masking-regime comparison
      Source: same CSV; dent=2 rows across all ablation modes
      Six regimes: full context, no geometry/stereo, ligand masked,
      ligand collapsed, metal+CN only, shuffled ligand (control)
      Grouped bar: per-donor Top-1 vs ligand exact match
      Orange bar = full-context exact match (reference); all others blue-grey
      Shuffled-ligand control confirms near-chance performance.
""")
(OUTDIR / "Fig_S7_notes.txt").write_text(notes)
print("Done: Fig_S7")
