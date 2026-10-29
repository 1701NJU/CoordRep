#!/usr/bin/env python3
"""
build_merged_fig5.py
====================
Generate figure-ready data for the merged Fig. 5:
  A. CoordRep-PathFinder workflow
  B. Full-CSD scope audit
  C. CN5 TBPY–SPY geometry ridge
  D. Family-level structural polymorphism
  E. Casebook mini-cards
  F. Quantitative output summary table

All numbers come from the full-CSD run (1,413,222 entries).
No CSD raw coordinates / CIF / XYZ exported.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

BASE = Path(__file__).parent.parent
SRC = BASE / "revision_results/csd_pathfinder_full"
OUT = BASE / "revision_results/figure_ready_merged_fig5"
os.makedirs(OUT, exist_ok=True)

LICENSE = "CSD-derived analysis results only; raw coordinates are not redistributed."

# ── Load source data ─────────────────────────────────────────────────

scan = json.load(open(SRC / "full_csd_scan_summary.json"))
metrics = json.load(open(SRC / "cn5_pathway_metrics.json"))
app_summary = json.load(open(SRC / "pathfinder_application_summary.json"))

TOTAL = scan["total_scanned"]          # 1,413,222
RETAINED = scan["retained"]            # 124,837
RETENTION = round(RETAINED / TOTAL * 100, 2)  # 8.83
N_BOUNDARY = scan["n_boundary"]        # 17,791
BOUNDARY_FRAC = round(N_BOUNDARY / RETAINED * 100, 1)  # 14.3
WATERFALL = scan["waterfall"]
REJECTIONS = scan["top_rejections"]

N_CN4_ATLAS = 47187
N_CN5_ATLAS = metrics["n_CN5"]         # 14,897
N_CN6_ATLAS = 41381
CN5_INTER_FRAC = round(metrics["fraction_delta_lt_1"] * 100, 1)  # 42.2
CN5_ENRICHMENT = round(metrics["ridge_enrichment_vs_random"], 1)  # 4.4

N_FAMILIES = app_summary["nontrivial_L3_families"]            # 6,358
N_MULTI_L1 = app_summary["families_with_multiple_L1"]         # 3,739
N_CROSS_BD = app_summary["families_crossing_shape_boundary"]  # 1,810

print("Source data loaded.")
print(f"  CSD: {TOTAL:,} scanned → {RETAINED:,} retained ({RETENTION}%)")
print(f"  CN5 atlas: {N_CN5_ATLAS:,}, intermediate: {CN5_INTER_FRAC}%")
print(f"  Families: {N_FAMILIES:,}, multi-L1: {N_MULTI_L1:,}, boundary-crossing: {N_CROSS_BD:,}")

# ══════════════════════════════════════════════════════════════════════
# Panel A: Workflow summary
# ══════════════════════════════════════════════════════════════════════

fig5A = {
    "input": "CSD structure",
    "record_layer": "CoordRep record",
    "derived_layers": ["CoordRep-ID", "validator", "CShM geometry fields"],
    "pathfinder_outputs": [
        "geometry atlas",
        "family trajectories",
        "boundary mining",
        "semantic audit"
    ],
    "caption_message": (
        "CSD structures are converted into auditable CoordRep records "
        "and then projected into geometry-pathway, identity-family, "
        "and curation outputs."
    ),
}
json.dump(fig5A, open(OUT / "fig5A_workflow_summary.json", "w"), indent=2)
print("\n[A] Workflow summary written.")


# ══════════════════════════════════════════════════════════════════════
# Panel B: Full-CSD scope audit
# ══════════════════════════════════════════════════════════════════════

# -- B.1: scope audit CSV --
with open(OUT / "fig5B_full_csd_scope_audit.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value", "label_for_figure"])
    w.writerow(["csd_release", "2024.3", "CSD release"])
    w.writerow(["total_entries_scanned", TOTAL,
                f"{TOTAL:,} CSD entries scanned"])
    w.writerow(["valid_coordrep_v1_records", RETAINED,
                f"{RETAINED:,} valid CoordRep v1 records"])
    w.writerow(["retention_rate_percent", RETENTION,
                f"{RETENTION}% retained"])

# -- B.2: scope audit JSON --
fig5B_json = {
    "csd_release": "2024.3",
    "total_entries_scanned": TOTAL,
    "valid_coordrep_v1_records": RETAINED,
    "retention_rate_percent": RETENTION,
    "waterfall_stages": {
        "has_3d": WATERFALL["has_3d"],
        "transition_metal": WATERFALL["transition_metal"],
        "mononuclear": WATERFALL["mononuclear"],
        "eta1": WATERFALL["eta1"],
        "cn_ok": WATERFALL["cn_ok"],
        "no_disorder": WATERFALL["no_disorder"],
        "no_polymer": WATERFALL["no_polymer"],
        "donor_ok": WATERFALL["donor_ok"],
        "valid_smiles": WATERFALL["valid_smiles"],
        "valid_coordrep": WATERFALL["valid_coordrep"],
    },
    "license_note": LICENSE,
}
json.dump(fig5B_json, open(OUT / "fig5B_full_csd_scope_audit.json", "w"), indent=2)

# -- B.3: outside-scope categories --
# Compute rejection counts from waterfall differences
total_rejected = TOTAL - RETAINED

# Readable label mapping
LABEL_MAP = {
    "no_TM":         "no transition metal",
    "multinuclear":  "multinuclear",
    "no_3d":         "no compatible 3D",
    "disorder":      "disorder / partial occupancy",
    "hapticity":     "unsupported ηn / haptic / non-η1",
    "raw_mol_fail":  "conversion failure",
    "polymeric":     "polymeric",
    "invalid_coordrep": "invalid CoordRep string",
}

# Take top 6 rejection reasons
sorted_rej = sorted(REJECTIONS.items(), key=lambda x: -x[1])

with open(OUT / "fig5B_outside_scope_categories.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["reason", "count", "percent_of_total_scanned",
                "percent_of_rejected", "label_for_figure"])
    for key, count in sorted_rej[:6]:
        label = LABEL_MAP.get(key, key)
        pct_total = round(count / TOTAL * 100, 2)
        pct_rej = round(count / total_rejected * 100, 2)
        w.writerow([key, count, pct_total, pct_rej, label])

print("[B] Scope audit files written.")


# ══════════════════════════════════════════════════════════════════════
# Panel C: CN5 TBPY–SPY geometry ridge
# ══════════════════════════════════════════════════════════════════════

# -- C.1: plot CSV from atlas --
atlas_path = SRC / "cn5_tbpy_spy_atlas.csv"
plot_rows = []
with open(atlas_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        s_tbpy = float(row["S_TBPY"])
        s_spy = float(row["S_SPY"])
        delta = abs(s_tbpy - s_spy)
        log_tbpy = math.log10(s_tbpy + 0.1)
        log_spy = math.log10(s_spy + 0.1)
        plot_rows.append({
            "refcode": row["refcode"],
            "metal": row["metal"],
            "oxidation_state": row["oxidation_state"],
            "d_count": row["d_count"],
            "donor_set": row["donor_set"],
            "denticity_pattern": row["denticity_pattern"],
            "S_TBPY": f"{s_tbpy:.4f}",
            "S_SPY": f"{s_spy:.4f}",
            "log_S_TBPY": f"{log_tbpy:.4f}",
            "log_S_SPY": f"{log_spy:.4f}",
            "delta_TBPY_SPY": f"{delta:.4f}",
            "is_intermediate_delta_lt_1": delta < 1.0,
            "L3_key_hash": row["L3_hash"],
        })

with open(OUT / "fig5C_cn5_tbpy_spy_plot.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=plot_rows[0].keys())
    w.writeheader()
    w.writerows(plot_rows)

# -- C.2: ridge annotations JSON --
fig5C_annotations = {
    "n_cn5_atlas_records": N_CN5_ATLAS,
    "intermediate_fraction_percent": CN5_INTER_FRAC,
    "criterion": "delta < 1.0",
    "enrichment_vs_random": CN5_ENRICHMENT,
    "ridge_metals": {
        "Cu": 1141,
        "Ru": 744,
        "Zn": 648,
        "Fe": 553,
        "V": 443,
    },
    "ridge_enrichment": {
        "monodentate_OR": 3.09,
        "C_N_P_donors_OR": 2.52,
        "Ir_Re_5d_OR": ">1.8",
    },
}
json.dump(fig5C_annotations, open(OUT / "fig5C_cn5_ridge_annotations.json", "w"), indent=2)

# -- C.3: ridge enrichment CSV --
enrichment_rows = [
    {"feature": "denticity_pattern", "category": "monodentate (1-1-1-1)",
     "count_in_ridge": 532, "odds_ratio": 3.09,
     "label_for_figure": "monodentate ligands OR=3.09"},
    {"feature": "donor_set", "category": "C,N,P",
     "count_in_ridge": 129, "odds_ratio": 2.52,
     "label_for_figure": "C,N,P donors OR=2.52"},
    {"feature": "metal", "category": "Ir",
     "count_in_ridge": 273, "odds_ratio": 2.46,
     "label_for_figure": "Ir OR=2.46"},
    {"feature": "metal", "category": "Ru",
     "count_in_ridge": 744, "odds_ratio": 2.32,
     "label_for_figure": "Ru OR=2.32"},
    {"feature": "metal_row", "category": "5d metals",
     "count_in_ridge": 1043, "odds_ratio": 1.72,
     "label_for_figure": "5d metals (row 6) OR=1.72"},
    {"feature": "metal", "category": "Re",
     "count_in_ridge": 257, "odds_ratio": 1.85,
     "label_for_figure": "Re OR=1.85"},
    {"feature": "metal", "category": "V",
     "count_in_ridge": 443, "odds_ratio": 1.67,
     "label_for_figure": "V OR=1.67"},
]
with open(OUT / "fig5C_cn5_ridge_enrichment.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=enrichment_rows[0].keys())
    w.writeheader()
    w.writerows(enrichment_rows)

print(f"[C] CN5 ridge files written ({len(plot_rows)} atlas records).")


# ══════════════════════════════════════════════════════════════════════
# Panel D: Family-level structural polymorphism
# ══════════════════════════════════════════════════════════════════════

# -- D.1: summary CSV --
with open(OUT / "fig5D_family_polymorphism_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["category", "count", "definition", "label_for_figure"])
    w.writerow(["nontrivial_L3_families", N_FAMILIES,
                "size >= 2 and unique L0 >= 2",
                f"{N_FAMILIES:,} nontrivial L3 families"])
    w.writerow(["multiple_L1_shape_families", N_MULTI_L1,
                "same L3 family contains multiple L1 ShapeIDs",
                f"{N_MULTI_L1:,} families with multiple L1 shapes"])
    w.writerow(["shape_boundary_crossing_families", N_CROSS_BD,
                "same L3 family crosses a shape-boundary criterion",
                f"{N_CROSS_BD:,} families crossing shape boundary"])

# -- D.2: top trajectory case JSON --
fig5D_case = {
    "case_name": "LINMOL",
    "metal": "Pd",
    "CN": 4,
    "trajectory_span": 10.10,
    "message": "Top geometry-trajectory span among full-CSD L3 families.",
    "n_records": 3,
    "n_refcodes": 3,
    "n_unique_L0": 3,
    "n_unique_L1": 3,
    "representative_refcodes": ["LINMOL", "LINMOL01", "LINMOL02"],
    "L3_key_hash": 463888590564,
    "crosses_shape_boundary": True,
}
json.dump(fig5D_case, open(OUT / "fig5D_top_trajectory_case.json", "w"), indent=2)

print("[D] Family polymorphism files written.")


# ══════════════════════════════════════════════════════════════════════
# Panel E: Casebook mini-cards
# ══════════════════════════════════════════════════════════════════════

cards = [
    {
        "case_id": "AFOSIA",
        "operation": "boundary geometry",
        "metal": "Cr",
        "CN": 6,
        "evidence_1": "Oh = 4.48",
        "evidence_2": "TPr = 4.49",
        "evidence_3": "Δ = 0.01",
        "key_result": "boundary",
        "figure_message": "top-2 shapes nearly tied",
    },
    {
        "case_id": "AGOTIA",
        "operation": "grammar repair",
        "metal": "Re",
        "CN": 7,
        "evidence_1": "missing bracket at pos. 81",
        "evidence_2": "invalid → valid",
        "evidence_3": "5/6 fields preserved",
        "key_result": "deterministic repair",
        "figure_message": "syntax by rules",
    },
    {
        "case_id": "CIJWUO",
        "operation": "stereo semantic check",
        "metal": "Cu",
        "CN": 5,
        "evidence_1": "trans → cis decoy",
        "evidence_2": "grammar-valid",
        "evidence_3": "margin = 0.10–0.68",
        "key_result": "Ranker 5/5 correct",
        "figure_message": "field-level consistency check",
    },
]

with open(OUT / "fig5E_casebook_minicards.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cards[0].keys())
    w.writeheader()
    w.writerows(cards)

cards_json = [dict(c, thumbnail_path="") for c in cards]
json.dump(cards_json, open(OUT / "fig5E_casebook_minicards.json", "w"), indent=2)

print("[E] Casebook mini-cards written.")


# ══════════════════════════════════════════════════════════════════════
# Panel F: Quantitative output summary table
# ══════════════════════════════════════════════════════════════════════

table_rows = [
    {
        "application": "Scope audit",
        "full_csd_output": f"{TOTAL:,} scanned; {RETAINED:,} valid",
        "chemical_meaning": "deterministic v1 domain",
        "representative_evidence": "Fig. 5B",
    },
    {
        "application": "CN5 geometry ridge",
        "full_csd_output": f"{N_CN5_ATLAS:,} records; {CN5_INTER_FRAC}% intermediate",
        "chemical_meaning": "static structures populate an interconversion coordinate",
        "representative_evidence": "Fig. 5C",
    },
    {
        "application": "Family trajectories",
        "full_csd_output": f"{N_FAMILIES:,} nontrivial L3 families",
        "chemical_meaning": "same connectivity spans multiple L0 states",
        "representative_evidence": "Fig. 5D",
    },
    {
        "application": "Shape polymorphism",
        "full_csd_output": f"{N_MULTI_L1:,} multi-L1 families; {N_CROSS_BD:,} crossing boundary",
        "chemical_meaning": "family-level geometry diversity",
        "representative_evidence": "Fig. 5D",
    },
    {
        "application": "Boundary mining",
        "full_csd_output": f"{N_BOUNDARY:,} records; {BOUNDARY_FRAC}%",
        "chemical_meaning": "boundary geometries are systematic",
        "representative_evidence": "Fig. 5E / AFOSIA",
    },
]

with open(OUT / "fig5F_quantitative_summary_table.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=table_rows[0].keys())
    w.writeheader()
    w.writerows(table_rows)

print("[F] Quantitative summary table written.")


# ══════════════════════════════════════════════════════════════════════
# All-caption numbers JSON
# ══════════════════════════════════════════════════════════════════════

caption = {
    "csd_release": "2024.3",
    "total_entries_scanned": TOTAL,
    "valid_coordrep_records": RETAINED,
    "retention_rate_percent": RETENTION,
    "cn4_atlas_records": N_CN4_ATLAS,
    "cn5_atlas_records": N_CN5_ATLAS,
    "cn6_atlas_records": N_CN6_ATLAS,
    "cn5_intermediate_fraction_percent": CN5_INTER_FRAC,
    "cn5_enrichment_vs_random": CN5_ENRICHMENT,
    "cn5_ridge_metals": {
        "Cu": 1141,
        "Ru": 744,
        "Zn": 648,
        "Fe": 553,
        "V": 443,
    },
    "cn5_monodentate_OR": 3.09,
    "cn5_CNP_donor_OR": 2.52,
    "nontrivial_L3_families": N_FAMILIES,
    "multiple_L1_shape_families": N_MULTI_L1,
    "shape_boundary_crossing_families": N_CROSS_BD,
    "top_trajectory_case": "Pd/CN4 LINMOL",
    "top_trajectory_span": 10.10,
    "boundary_records": N_BOUNDARY,
    "boundary_fraction_percent": BOUNDARY_FRAC,
    "casebook": {
        "AFOSIA": "Oh 4.48 vs TPr 4.49; delta 0.01",
        "AGOTIA": "missing bracket; invalid to valid; 5/6 fields preserved",
        "CIJWUO": "trans to cis decoy; 5/5 correct; margin 0.10–0.68",
    },
}
json.dump(caption, open(OUT / "fig5_all_caption_numbers.json", "w"), indent=2)

print("[✓] All-caption numbers written.")


# ══════════════════════════════════════════════════════════════════════
# README.md
# ══════════════════════════════════════════════════════════════════════

readme = f"""# Figure 5 — Figure-Ready Data (Merged)

**Figure 5.** Full-CSD CoordRep-PathFinder maps coordination-geometry continua
and record-level structural diversity.

## Panels

| Panel | Title | Key file |
|-------|-------|----------|
| A | CoordRep-PathFinder workflow | `fig5A_workflow_summary.json` |
| B | Full-CSD scope audit | `fig5B_full_csd_scope_audit.csv` |
| C | CN5 TBPY–SPY geometry ridge | `fig5C_cn5_tbpy_spy_plot.csv` |
| D | Family-level structural polymorphism | `fig5D_family_polymorphism_summary.csv` |
| E | Casebook mini-cards | `fig5E_casebook_minicards.csv` |
| F | Quantitative output summary table | `fig5F_quantitative_summary_table.csv` |

## Data provenance

1. **All Fig. 5 numbers come from the full-CSD run** (1,413,222 entries scanned).
   The old 200K sample is no longer used.
2. **Fig. 5B** rejection categories come from the full-CSD filtering waterfall
   (`full_csd_scan_summary.json → waterfall`).
3. **Fig. 5C** CN5 atlas records: **{N_CN5_ATLAS:,}**.
4. **Fig. 5D** family counts:
   - {N_FAMILIES:,} nontrivial L3 families
   - {N_MULTI_L1:,} multiple L1 families
   - {N_CROSS_BD:,} crossing boundary families
5. **Boundary records**: {N_BOUNDARY:,} ({BOUNDARY_FRAC}%).
6. **No CSD raw coordinates / CIF / XYZ exported.**
7. Output files contain only refcodes, hashed keys, aggregate statistics,
   CShM values, and field summaries.

## Key numbers (for caption cross-check)

| Metric | Value |
|--------|-------|
| CSD release | 2024.3 |
| Total entries scanned | {TOTAL:,} |
| Valid CoordRep v1 records | {RETAINED:,} |
| Retention rate | {RETENTION}% |
| CN4 atlas | {N_CN4_ATLAS:,} |
| CN5 atlas | {N_CN5_ATLAS:,} |
| CN6 atlas | {N_CN6_ATLAS:,} |
| CN5 intermediate (delta<1) | {CN5_INTER_FRAC}% |
| CN5 enrichment vs random | {CN5_ENRICHMENT}× |
| Nontrivial L3 families | {N_FAMILIES:,} |
| Multi-L1 families | {N_MULTI_L1:,} |
| Boundary-crossing families | {N_CROSS_BD:,} |
| Boundary records | {N_BOUNDARY:,} ({BOUNDARY_FRAC}%) |
| Top trajectory | Pd/CN4 LINMOL, span = 10.10 |

## License

{LICENSE}

## Files

```
figure_ready_merged_fig5/
├── README.md
├── fig5A_workflow_summary.json
├── fig5B_full_csd_scope_audit.csv
├── fig5B_full_csd_scope_audit.json
├── fig5B_outside_scope_categories.csv
├── fig5C_cn5_tbpy_spy_plot.csv
├── fig5C_cn5_ridge_annotations.json
├── fig5C_cn5_ridge_enrichment.csv
├── fig5D_family_polymorphism_summary.csv
├── fig5D_top_trajectory_case.json
├── fig5E_casebook_minicards.csv
├── fig5E_casebook_minicards.json
├── fig5F_quantitative_summary_table.csv
├── fig5_all_caption_numbers.json
├── fig5B_scope_audit_preview.png          (optional)
├── fig5C_cn5_tbpy_spy_hexbin_preview.png  (optional)
├── fig5D_family_polymorphism_bar_preview.png (optional)
└── fig5F_summary_table_preview.png        (optional)
```
"""
with open(OUT / "README.md", "w") as f:
    f.write(readme)

print("[✓] README.md written.")


# ══════════════════════════════════════════════════════════════════════
# Quick-look preview plots
# ══════════════════════════════════════════════════════════════════════

ORANGE = "#E87D2F"
GRAY_DARK = "#333333"
GRAY_MED = "#888888"
GRAY_LIGHT = "#CCCCCC"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.edgecolor": GRAY_DARK,
    "axes.linewidth": 0.8,
    "xtick.color": GRAY_DARK,
    "ytick.color": GRAY_DARK,
})


# ── Fig 5B: Scope audit waterfall ────────────────────────────────────

def plot_5B():
    stages = [
        ("CSD entries", TOTAL),
        ("3D structure", WATERFALL["has_3d"]),
        ("transition metal", WATERFALL["transition_metal"]),
        ("mononuclear", WATERFALL["mononuclear"]),
        ("η1 only", WATERFALL["eta1"]),
        ("CN 2–6", WATERFALL["cn_ok"]),
        ("no disorder", WATERFALL["no_disorder"]),
        ("no polymer", WATERFALL["no_polymer"]),
        ("valid CoordRep", WATERFALL["valid_coordrep"]),
    ]
    labels = [s[0] for s in stages]
    values = [s[1] for s in stages]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(range(len(labels)), values, color=GRAY_LIGHT, edgecolor=GRAY_MED)
    bars[-1].set_color(ORANGE)
    bars[-1].set_edgecolor(ORANGE)

    for i, (lab, val) in enumerate(stages):
        ax.text(val + TOTAL * 0.01, i, f"{val:,}", va="center", fontsize=8,
                color=GRAY_DARK)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Number of entries")
    ax.set_title("Fig. 5B — Full-CSD Scope Audit", fontsize=11, fontweight="bold",
                 color=GRAY_DARK)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}M"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig5B_scope_audit_preview.png", dpi=200, facecolor="white")
    plt.close(fig)
    print("[plot] fig5B_scope_audit_preview.png")


# ── Fig 5C: CN5 hexbin ──────────────────────────────────────────────

def plot_5C():
    xs, ys, is_inter = [], [], []
    for row in plot_rows:
        xs.append(float(row["log_S_TBPY"]))
        ys.append(float(row["log_S_SPY"]))
        is_inter.append(row["is_intermediate_delta_lt_1"])

    xs = np.array(xs)
    ys = np.array(ys)
    inter = np.array(is_inter)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    hb = ax.hexbin(xs, ys, gridsize=40, cmap="Greys", mincnt=1)
    cb = fig.colorbar(hb, ax=ax, shrink=0.7, label="Count")

    # Overlay ridge zone
    mn, mx = min(xs.min(), ys.min()) - 0.1, max(xs.max(), ys.max()) + 0.1
    ax.plot([mn, mx], [mn, mx], "--", color=ORANGE, lw=1.5, alpha=0.8,
            label="S_TBPY = S_SPY")

    ax.set_xlabel("log₁₀(S_TBPY + 0.1)")
    ax.set_ylabel("log₁₀(S_SPY + 0.1)")
    ax.set_title("Fig. 5C — CN5 TBPY–SPY Geometry Ridge", fontsize=11,
                 fontweight="bold", color=GRAY_DARK)

    # Annotation
    n_inter = int(inter.sum())
    ax.text(0.03, 0.97,
            f"n = {N_CN5_ATLAS:,}\n"
            f"Δ < 1.0: {CN5_INTER_FRAC}% ({n_inter:,})\n"
            f"enrichment: {CN5_ENRICHMENT}×",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(facecolor="white", edgecolor=GRAY_MED, alpha=0.9))
    ax.legend(loc="lower right", fontsize=8)
    ax.set_aspect("equal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig5C_cn5_tbpy_spy_hexbin_preview.png", dpi=200,
                facecolor="white")
    plt.close(fig)
    print("[plot] fig5C_cn5_tbpy_spy_hexbin_preview.png")


# ── Fig 5D: Family polymorphism bars ─────────────────────────────────

def plot_5D():
    cats = [
        f"Nontrivial L3\nfamilies",
        f"Multiple L1\nshapes",
        f"Crossing\nboundary",
    ]
    vals = [N_FAMILIES, N_MULTI_L1, N_CROSS_BD]

    fig, ax = plt.subplots(figsize=(5, 3.5))
    bars = ax.bar(range(3), vals, color=[GRAY_MED, GRAY_MED, ORANGE],
                  edgecolor=GRAY_DARK, width=0.6)

    for i, v in enumerate(vals):
        ax.text(i, v + 100, f"{v:,}", ha="center", fontsize=10,
                fontweight="bold", color=GRAY_DARK)

    ax.set_xticks(range(3))
    ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("Number of families")
    ax.set_title("Fig. 5D — Family-Level Structural Polymorphism",
                 fontsize=11, fontweight="bold", color=GRAY_DARK)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig5D_family_polymorphism_bar_preview.png", dpi=200,
                facecolor="white")
    plt.close(fig)
    print("[plot] fig5D_family_polymorphism_bar_preview.png")


# ── Fig 5F: Summary table ───────────────────────────────────────────

def plot_5F():
    col_labels = ["Application", "Full-CSD output", "Chemical meaning", "Evidence"]
    cell_text = []
    for r in table_rows:
        cell_text.append([
            r["application"],
            r["full_csd_output"],
            r["chemical_meaning"],
            r["representative_evidence"],
        ])

    fig, ax = plt.subplots(figsize=(10, 2.5))
    ax.axis("off")
    ax.set_title("Fig. 5F — Quantitative Output Summary", fontsize=11,
                 fontweight="bold", color=GRAY_DARK, pad=12)

    tbl = ax.table(cellText=cell_text, colLabels=col_labels,
                   loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.6)

    # Style header
    for j in range(len(col_labels)):
        cell = tbl[0, j]
        cell.set_facecolor(GRAY_DARK)
        cell.set_text_props(color="white", fontweight="bold")

    # Alternate row shading
    for i in range(1, len(cell_text) + 1):
        for j in range(len(col_labels)):
            cell = tbl[i, j]
            if i % 2 == 0:
                cell.set_facecolor("#F5F5F5")
            else:
                cell.set_facecolor("white")

    fig.tight_layout()
    fig.savefig(OUT / "fig5F_summary_table_preview.png", dpi=200,
                facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("[plot] fig5F_summary_table_preview.png")


plot_5B()
plot_5C()
plot_5D()
plot_5F()


# ══════════════════════════════════════════════════════════════════════
# Final check
# ══════════════════════════════════════════════════════════════════════

required = [
    "fig5A_workflow_summary.json",
    "fig5B_full_csd_scope_audit.csv",
    "fig5B_full_csd_scope_audit.json",
    "fig5B_outside_scope_categories.csv",
    "fig5C_cn5_tbpy_spy_plot.csv",
    "fig5C_cn5_ridge_annotations.json",
    "fig5C_cn5_ridge_enrichment.csv",
    "fig5D_family_polymorphism_summary.csv",
    "fig5D_top_trajectory_case.json",
    "fig5E_casebook_minicards.csv",
    "fig5E_casebook_minicards.json",
    "fig5F_quantitative_summary_table.csv",
    "fig5_all_caption_numbers.json",
    "README.md",
]

print(f"\n{'='*60}")
print("Output verification")
print(f"{'='*60}")
ok = 0
for fn in required:
    p = OUT / fn
    if p.exists():
        sz = os.path.getsize(p)
        print(f"  ✓ {fn:50s} {sz:>8,} B")
        ok += 1
    else:
        print(f"  ✗ {fn:50s} MISSING")

# Optional previews
for fn in ["fig5B_scope_audit_preview.png",
           "fig5C_cn5_tbpy_spy_hexbin_preview.png",
           "fig5D_family_polymorphism_bar_preview.png",
           "fig5F_summary_table_preview.png"]:
    p = OUT / fn
    if p.exists():
        sz = os.path.getsize(p)
        print(f"  ✓ {fn:50s} {sz:>8,} B  (preview)")

print(f"\nRequired: {ok}/{len(required)} files present.")
print(f"{'='*60}")
