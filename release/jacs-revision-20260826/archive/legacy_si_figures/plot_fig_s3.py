#!/usr/bin/env python3
"""
Supplementary Figure S3 — CoordRep-ID robustness across identity layers.

Input files:
  revision_results/identity_robustness/rigid_body_invariance_test.csv
  revision_results/identity_robustness/perturbation_sweep.csv
  revision_results/identity_robustness/csd_refcode_family_validation.csv
  revision_results/identity_robustness/bond_scaling_sweep.csv
  revision_results/identity_robustness/family_case_study.csv

Output files:
  revision_results/si_figures/Fig_S3.{pdf,svg,png}
  revision_results/si_figures/Fig_S3_notes.txt

Panels:
  A — Deterministic invariance tests (bar chart, L0 retention = 100%)
  B — Coordinate perturbation sweep (line plot, L0–L3 vs σ, log10 x)
  C — CSD refcode-family validation (empirical, bar chart)
  D — Bond-scaling robustness (line plot, L0–L3 vs scale factor)
"""
import sys, csv, textwrap
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "si_figures"))

from jacs_style import (apply_jacs_style, save_fig, label_panel, clean_axes,
                         light_grid, BLACK, DARK_GREY, MID_GREY, BLUE_GREY,
                         ACCENT, PALE_GREY, FAINT_GREY, LAYER_COLORS)
import matplotlib.pyplot as plt

apply_jacs_style()

ROBDIR = ROOT / "revision_results" / "identity_robustness"
OUTDIR = ROOT / "revision_results" / "si_figures"

# ── Load data ────────────────────────────────────────
def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

inv       = load_csv(ROBDIR / "rigid_body_invariance_test.csv")
sweep     = load_csv(ROBDIR / "perturbation_sweep.csv")
fam_val   = load_csv(ROBDIR / "csd_refcode_family_validation.csv")
bscale    = load_csv(ROBDIR / "bond_scaling_sweep.csv")
case_rows = load_csv(ROBDIR / "family_case_study.csv")

# ── Figure layout ────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.5))
fig.subplots_adjust(hspace=0.45, wspace=0.35,
                    left=0.10, right=0.97, top=0.95, bottom=0.08)
ax_a, ax_b, ax_c, ax_d = axes.flat

# ═══════════════════════════════════════════════════════
# Panel A — Deterministic invariance tests
# ═══════════════════════════════════════════════════════
label_panel(ax_a, "A")
clean_axes(ax_a)
light_grid(ax_a)

tests = [r["test"].replace("_", "\n") for r in inv]
vals = [float(r["L0_retention"]) * 100 for r in inv]
x = np.arange(len(tests))
bars = ax_a.bar(x, vals, width=0.55, color=DARK_GREY, edgecolor=BLACK, linewidth=0.4)
ax_a.set_xticks(x)
ax_a.set_xticklabels(tests, fontsize=6)
ax_a.set_ylabel("L0 retention (%)")
ax_a.set_ylim(98, 101)
ax_a.set_title("Deterministic invariance", fontsize=8)
for b, v in zip(bars, vals):
    ax_a.text(b.get_x() + b.get_width()/2, v + 0.2, f"{v:.0f}%",
              ha="center", va="bottom", fontsize=6, color=DARK_GREY)

# ═══════════════════════════════════════════════════════
# Panel B — Coordinate perturbation sweep (log10 x-axis)
# ═══════════════════════════════════════════════════════
label_panel(ax_b, "B")
clean_axes(ax_b)
light_grid(ax_b)

# Skip σ=0 for log scale
sweep_pos = [r for r in sweep if float(r["sigma_angstrom"]) > 0]
sigmas = [float(r["sigma_angstrom"]) for r in sweep_pos]
for level in ["L0", "L1", "L2", "L3"]:
    vals = [float(r[f"{level}_retention"]) * 100 for r in sweep_pos]
    ax_b.plot(sigmas, vals, marker="o", markersize=3.5,
              color=LAYER_COLORS[level], label=level, linewidth=1.3)

ax_b.set_xlabel("Perturbation σ (Å)")
ax_b.set_ylabel("Retention (%)")
ax_b.set_xscale("log")
ax_b.set_ylim(-2, 105)
ax_b.legend(loc="lower left", fontsize=6.5, ncol=2,
            columnspacing=0.8, handletextpad=0.4)
ax_b.set_title("Perturbation robustness", fontsize=8)

# ═══════════════════════════════════════════════════════
# Panel C — CSD refcode-family validation (empirical)
# ═══════════════════════════════════════════════════════
label_panel(ax_c, "C")
clean_axes(ax_c)
light_grid(ax_c)

levels_c = [r["level"] for r in fam_val]
match_rates = [float(r["within_family_match_rate"]) * 100 for r in fam_val]
colors_c = [LAYER_COLORS.get(l, MID_GREY) for l in levels_c]
n_fam = int(fam_val[0]["n_families_evaluated"])

x = np.arange(len(levels_c))
bars_c = ax_c.bar(x, match_rates, width=0.55, color=colors_c,
                   edgecolor=BLACK, linewidth=0.4)
ax_c.set_xticks(x)
ax_c.set_xticklabels(levels_c)
ax_c.set_ylabel("Within-family\nmatch rate (%)")
ax_c.set_ylim(0, 110)
ax_c.set_title("CSD refcode-family validation", fontsize=8)
ax_c.text(0.02, 0.95, f"n = {n_fam:,} CSD refcode families",
          transform=ax_c.transAxes, ha="left", va="top",
          fontsize=6, color=MID_GREY)
for b, v in zip(bars_c, match_rates):
    ax_c.text(b.get_x() + b.get_width()/2, v + 1.5, f"{v:.1f}%",
              ha="center", va="bottom", fontsize=6, color=DARK_GREY)

# Case study annotation — positioned below the bars
if case_rows:
    cs = case_rows[0]
    cs_text = (f"e.g. {cs['family']}: {cs['n_members']} members, "
               f"{cs['unique_L0']} L0 → {cs['unique_L3']} L3")
    ax_c.text(0.50, -0.18, cs_text, transform=ax_c.transAxes,
              ha="center", va="top", fontsize=5.5, color=ACCENT,
              fontstyle="italic")

# ═══════════════════════════════════════════════════════
# Panel D — Bond-scaling robustness
# ═══════════════════════════════════════════════════════
label_panel(ax_d, "D")
clean_axes(ax_d)
light_grid(ax_d)

scale_f = [float(r["scale_factor"]) for r in bscale]
for level, key in [("L0", "L0_retention"), ("L1", "L1_retention"),
                    ("L3", "L3_retention")]:
    vals = [float(r[key]) * 100 for r in bscale]
    ax_d.plot(scale_f, vals, marker="o", markersize=3.5,
              color=LAYER_COLORS[level], label=level, linewidth=1.3)

ax_d.axvline(1.0, ls=":", lw=0.6, color=PALE_GREY, zorder=0)
ax_d.set_xlabel("Bond scaling factor")
ax_d.set_ylabel("Retention (%)")
ax_d.set_ylim(90, 101)
ax_d.set_title("Bond-scaling robustness", fontsize=8)
ax_d.legend(loc="lower center", fontsize=6.5, ncol=3,
            columnspacing=0.8, handletextpad=0.4)

# ── Save ─────────────────────────────────────────────
save_fig(fig, "Fig_S3", OUTDIR)
plt.close(fig)

# ── Notes ────────────────────────────────────────────
cs = case_rows[0] if case_rows else {}
notes = textwrap.dedent(f"""\
    Fig_S3 — CoordRep-ID robustness across identity layers

    Panel A: Deterministic invariance tests
      Source: identity_robustness/rigid_body_invariance_test.csv
      x = test type; y = L0 retention (%)
      All nuisance transformations preserve L0, confirming deterministic
      canonicalization under the adopted rules.

    Panel B: Coordinate perturbation sweep
      Source: identity_robustness/perturbation_sweep.csv
      x = perturbation σ (Å), log10 scale; y = identity retention (%)
      Lines: L0 (dark grey), L1 (mid grey), L3 (accent orange)
      L0 is most sensitive because it preserves the exact geometry-resolved
      state, whereas L1–L3 progressively abstract continuous geometry.

    Panel C: CSD refcode-family validation
      Source: identity_robustness/csd_refcode_family_validation.csv
      Families defined by CSD refcode families (independent of CoordRep).
      x = identity level; y = within-family match rate (%)
      n = {n_fam:,} CSD refcode families with ≥2 members.
      L0 = 11.4%, L1 = 54.3%, L3 = 93.7%
      Demonstrates that coarser identity layers are more stable across
      independent experimental redeterminations.
      Case study: {cs.get('family','')} ({cs.get('n_members','')} members,
      {cs.get('unique_L0','')} unique L0 → {cs.get('unique_L3','')} shared L3)

    Panel D: Bond-scaling robustness
      Source: identity_robustness/bond_scaling_sweep.csv
      x = uniform bond scaling factor; y = retention (%)
      Lines: L0, L1, L3
      CShM is approximately scale-invariant, so higher-level keys are robust
      to uniform bond-length changes; L0 is affected by numerical CShM shifts.
""")
(OUTDIR / "Fig_S3_notes.txt").write_text(notes)
print("Done: Fig_S3")
