#!/usr/bin/env python3
"""
Generate CSD scope coverage waterfall for SI.

Responds to the reviewer question: "What fraction of real coordination
complexes are excluded by current filtering criteria?"

Instead of reporting only 124,837 / 1,413,222 = 8.83%, this script
produces a multi-denominator waterfall with hierarchical exclusion
categories and denominator-corrected metrics.

Input:  CSD 2024.3 full-scan filtering logs
Output: revision_results/csd_scope_waterfall_revision/
"""

import csv
import json
import textwrap
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "revision_results" / "csd_scope_waterfall_revision"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ════════════════════════════════════════════════════════════════════
# Source data — from full_csd_scan_summary.json (CSD 2024.3)
# ════════════════════════════════════════════════════════════════════
# These are the cumulative waterfall counts from the sequential filter
# pipeline.  Each stage is the number of entries surviving UP TO that
# filter.

WATERFALL_CUMULATIVE = OrderedDict([
    ("total_scanned",       1_413_222),
    ("has_3d",              1_335_141),
    ("transition_metal",      615_498),
    ("mononuclear",           269_030),
    ("eta1",                  217_163),
    ("cn_ok",                 209_826),
    ("no_disorder",           145_728),
    ("no_polymer",            144_271),
    ("donor_ok",              126_197),
    ("valid_smiles",          126_197),
    ("valid_coordrep",        124_837),
])

# top_rejections from the same scan (primary reason per entry)
TOP_REJECTIONS = {
    "no_TM":            719_643,
    "multinuclear":     346_468,
    "no_3d":             78_081,
    "disorder":          64_098,
    "hapticity":         51_867,
    "raw_mol_fail":      18_074,
    "cn_7":               3_338,
    "cn_8":               1_823,
    "polymeric":          1_457,
    "invalid_coordrep":   1_360,
    "cn_9":                 536,
    "cn_0":                 510,
    "cn_1":                 412,
    "cn_10":                346,
    "cn_12":                196,
}

# ════════════════════════════════════════════════════════════════════
# Derived counts — mutually exclusive waterfall layers
# ════════════════════════════════════════════════════════════════════
# Priority order (highest to lowest):
#   1. no transition metal at all
#   2. no compatible 3D / atom-resolved structure
#   3. → transition_metal_candidate_entries (checkpoint)
#   4. multinuclear / extended coordination
#   5. haptic / π (η > 1)
#   6. disorder / partial occupancy
#   7. incompatible or ambiguous bonding graph (CN out of range,
#      polymer, raw-mol failure, donor-fail)
#   8. outside v1 grammar — other (residual)
#   9. → intended_v1_candidate_entries (checkpoint)
#  10. valid CoordRep v1 records
#  11. in-domain conversion failures

W = WATERFALL_CUMULATIVE
total = W["total_scanned"]

# ── Mutually exclusive exclusion categories ──────────────────────
# We apply the same sequential priority as the pipeline.
# Entries removed at an earlier stage are NOT counted again later.

# no_transition_metal: entries that had 3D but no TM, PLUS entries
# with no 3D (since they also lack TM).  But to keep things clean:
# "no compatible 3D" = total - has_3d
# "no TM" = has_3d - transition_metal (these had 3D but no TM)
no_3d = total - W["has_3d"]                               # 78,081
no_tm = W["has_3d"] - W["transition_metal"]                # 719,643

# TM candidates (checkpoint)
tm_candidates = W["transition_metal"]                      # 615,498

# multinuclear: TM entries that are not mononuclear
multinuclear = W["transition_metal"] - W["mononuclear"]    # 346,468

# haptic / π: mononuclear but η > 1
haptic_pi = W["mononuclear"] - W["eta1"]                   # 51,867

# CN out of range (mononuclear, η1, but CN outside 2–6)
cn_out = W["eta1"] - W["cn_ok"]                            # 7,337

# disorder / partial occupancy
disorder = W["cn_ok"] - W["no_disorder"]                   # 64,098

# polymer
polymer = W["no_disorder"] - W["no_polymer"]               # 1,457

# donor / bonding / SMILES failures
donor_fail = W["no_polymer"] - W["valid_smiles"]           # 18,074

# Aggregate: incompatible or ambiguous bonding graph
incompatible_bonding = cn_out + polymer + donor_fail       # 26,868

# Outside current v1 grammar — other residual
# (any entries not yet categorized)
# All entries accounted for by the above + disorder + valid, so other = 0
# Let's verify:
accounted = no_3d + no_tm + multinuclear + haptic_pi + incompatible_bonding + disorder
intended_v1_candidates = total - accounted
# intended_v1_candidates should equal valid_smiles (126,197) ideally,
# but let's compute precisely
# Actually: intended_v1 = entries surviving all exclusions except
# final coordrep validation
# From waterfall: valid_smiles = 126,197 (same as donor_ok)
# valid_coordrep = 124,837
# So in-domain conversion failures = 126,197 - 124,837 = 1,360

in_domain_failures = W["valid_smiles"] - W["valid_coordrep"]  # 1,360
valid_records = W["valid_coordrep"]                            # 124,837

# outside_v1_other: anything not yet categorized
outside_v1_other = total - (no_3d + no_tm + multinuclear + haptic_pi
                            + incompatible_bonding + disorder
                            + in_domain_failures + valid_records)

# Verify balance
assert no_3d + no_tm + multinuclear + haptic_pi + incompatible_bonding \
    + disorder + outside_v1_other + in_domain_failures + valid_records == total, \
    "Waterfall does not sum to total!"

# ── Intermediate checkpoints ──
# coordination_candidate = entries with TM, mononuclear, η1
coordination_candidates = W["eta1"]                        # 217,163

# intended_v1_candidate = entries surviving all exclusions except
# final validation
intended_v1 = valid_records + in_domain_failures           # 126,197


# ════════════════════════════════════════════════════════════════════
# A. Hierarchical waterfall CSV
# ════════════════════════════════════════════════════════════════════

def write_waterfall_csv():
    rows = [
        ("total_CSD_entries",                           total,
         "",    "Total CSD 2024.3 entries"),
        ("no_compatible_3D_or_atom_resolved_structure",  no_3d,
         "excluded", "No 3D coordinates or atom-resolved structure"),
        ("no_transition_metal",                          no_tm,
         "excluded", "3D present but no transition metal"),
        ("transition_metal_candidate_entries",            tm_candidates,
         "checkpoint", "Entries containing ≥1 transition metal with 3D"),
        ("multinuclear_or_extended_coordination",        multinuclear,
         "excluded", "≥2 metal centers (dinuclear, cluster, MOF, etc.)"),
        ("haptic_or_pi_coordination_eta_greater_than_1", haptic_pi,
         "excluded", "Mononuclear but with η>1 haptic/π coordination"),
        ("disorder_or_partial_occupancy",                disorder,
         "excluded", "Crystallographic disorder or partial occupancy"),
        ("incompatible_or_ambiguous_bonding_graph",      incompatible_bonding,
         "excluded", "CN outside 2–6, polymeric, or donor/SMILES failure"),
        ("outside_current_v1_grammar_other",             outside_v1_other,
         "excluded", "Other residual exclusions"),
        ("intended_v1_candidate_entries",                 intended_v1,
         "checkpoint", "Entries within intended CoordRep v1 domain"),
        ("in_domain_conversion_failures",                in_domain_failures,
         "excluded", "In-domain but failed final CoordRep validation"),
        ("valid_CoordRep_v1_records",                    valid_records,
         "retained", "Successfully converted to valid CoordRep v1"),
    ]

    with open(OUTDIR / "csd_scope_waterfall.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "n_entries", "type", "description"])
        for row in rows:
            w.writerow(row)

    print(f"  csd_scope_waterfall.csv  ({len(rows)} rows)")
    return rows


# ════════════════════════════════════════════════════════════════════
# B. Summary JSON with denominator-corrected metrics
# ════════════════════════════════════════════════════════════════════

def write_summary_json():
    summary = OrderedDict([
        ("csd_release",                                  "2024.3"),
        ("total_CSD_entries",                            total),
        ("valid_CoordRep_v1_records",                    valid_records),
        ("whole_CSD_fraction_percent",                   round(valid_records / total * 100, 2)),

        ("transition_metal_candidate_entries",            tm_candidates),
        ("coordination_candidate_entries",                coordination_candidates),
        ("intended_v1_candidate_entries",                 intended_v1),

        ("valid_fraction_among_TM_candidates_percent",
         round(valid_records / tm_candidates * 100, 2)),
        ("valid_fraction_among_coordination_candidates_percent",
         round(valid_records / coordination_candidates * 100, 2)),
        ("valid_fraction_among_intended_v1_candidates_percent",
         round(valid_records / intended_v1 * 100, 2)),

        ("excluded_fraction_among_coordination_candidates_percent",
         round((coordination_candidates - valid_records) / coordination_candidates * 100, 2)),
        ("excluded_fraction_due_to_multinuclear_percent",
         round(multinuclear / tm_candidates * 100, 2)),
        ("excluded_fraction_due_to_haptic_or_pi_percent",
         round(haptic_pi / tm_candidates * 100, 2)),
        ("excluded_fraction_due_to_disorder_percent",
         round(disorder / tm_candidates * 100, 2)),
        ("excluded_fraction_due_to_incompatible_bonding_percent",
         round(incompatible_bonding / tm_candidates * 100, 2)),
        ("excluded_fraction_due_to_incompatible_3D_percent",
         round(no_3d / total * 100, 2)),

        ("exclusion_priority_note",
         "Each entry is classified by its first (highest-priority) "
         "exclusion reason.  Priority: no-3D > no-TM > multinuclear > "
         "haptic/π > CN-out-of-range > disorder > polymer > "
         "donor/SMILES-fail > other > in-domain-validation-fail."),
        ("license_note",
         "CSD-derived aggregate counts only; no raw coordinates are "
         "redistributed."),
    ])

    with open(OUTDIR / "csd_scope_waterfall_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("  csd_scope_waterfall_summary.json")
    return summary


# ════════════════════════════════════════════════════════════════════
# C. Exclusion-by-denominator table
# ════════════════════════════════════════════════════════════════════

def write_exclusion_by_denominator():
    categories = [
        ("total_CSD_entries",                           total,
         "Total entries in CSD 2024.3"),
        ("no_compatible_3D",                            no_3d,
         "Lack atom-resolved 3D coordinates"),
        ("no_transition_metal",                         no_tm,
         "Organic or main-group only"),
        ("transition_metal_candidates",                  tm_candidates,
         "Contain ≥1 TM with 3D (denominator 2)"),
        ("multinuclear_or_extended",                    multinuclear,
         "Di-/polynuclear, clusters, MOFs"),
        ("haptic_or_pi",                                haptic_pi,
         "Mononuclear with η>1 haptic/π sites"),
        ("disorder_or_partial_occupancy",               disorder,
         "Crystallographic disorder prevents unambiguous encoding"),
        ("incompatible_bonding_graph",                  incompatible_bonding,
         "CN outside 2–6, polymeric, donor/SMILES issues"),
        ("outside_v1_grammar_other",                    outside_v1_other,
         "Residual exclusions"),
        ("intended_v1_candidates",                       intended_v1,
         "Within CoordRep v1 intended scope (denominator 3)"),
        ("in_domain_conversion_failures",               in_domain_failures,
         "In-domain entries failing final serialization"),
        ("valid_CoordRep_v1_records",                   valid_records,
         "Final valid records"),
    ]

    def pct(n, d):
        return round(n / d * 100, 2) if d > 0 else 0.0

    fields = ["category", "n_entries", "percent_of_total_CSD",
              "percent_of_TM_candidates", "percent_of_coordination_candidates",
              "percent_of_v1_candidates", "interpretation"]

    with open(OUTDIR / "csd_scope_exclusion_by_denominator.csv",
              "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, n, interp in categories:
            w.writerow({
                "category": name,
                "n_entries": n,
                "percent_of_total_CSD": pct(n, total),
                "percent_of_TM_candidates": pct(n, tm_candidates),
                "percent_of_coordination_candidates": pct(n, coordination_candidates),
                "percent_of_v1_candidates": pct(n, intended_v1),
                "interpretation": interp,
            })

    print("  csd_scope_exclusion_by_denominator.csv")


# ════════════════════════════════════════════════════════════════════
# D. Figure 5B updated scope waterfall data
# ════════════════════════════════════════════════════════════════════

def write_fig5b_data():
    rows = [
        {
            "layer": "All CSD entries",
            "n_entries": total,
            "label": f"{total:,} entries (CSD 2024.3)",
            "fraction_of_previous": 1.0,
            "fraction_of_total": 1.0,
            "bar_color": "grey",
        },
        {
            "layer": "Transition-metal candidates",
            "n_entries": tm_candidates,
            "label": f"{tm_candidates:,} TM entries",
            "fraction_of_previous": round(tm_candidates / total, 4),
            "fraction_of_total": round(tm_candidates / total, 4),
            "bar_color": "blue_grey",
        },
        {
            "layer": "Mononuclear η1 coordination candidates",
            "n_entries": coordination_candidates,
            "label": f"{coordination_candidates:,} coordination candidates",
            "fraction_of_previous": round(coordination_candidates / tm_candidates, 4),
            "fraction_of_total": round(coordination_candidates / total, 4),
            "bar_color": "blue_grey",
        },
        {
            "layer": "Intended CoordRep v1 domain",
            "n_entries": intended_v1,
            "label": f"{intended_v1:,} v1 candidates",
            "fraction_of_previous": round(intended_v1 / coordination_candidates, 4),
            "fraction_of_total": round(intended_v1 / total, 4),
            "bar_color": "blue_grey",
        },
        {
            "layer": "Valid CoordRep v1 records",
            "n_entries": valid_records,
            "label": f"{valid_records:,} valid records",
            "fraction_of_previous": round(valid_records / intended_v1, 4),
            "fraction_of_total": round(valid_records / total, 4),
            "bar_color": "orange",
        },
    ]

    fields = list(rows[0].keys())
    with open(OUTDIR / "fig5B_updated_scope_waterfall.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print("  fig5B_updated_scope_waterfall.csv")


# ════════════════════════════════════════════════════════════════════
# E. SI markdown + README
# ════════════════════════════════════════════════════════════════════

def write_si_markdown():
    pct_tm = round(valid_records / tm_candidates * 100, 1)
    pct_coord = round(valid_records / coordination_candidates * 100, 1)
    pct_v1 = round(valid_records / intended_v1 * 100, 1)
    pct_multi = round(multinuclear / tm_candidates * 100, 1)
    pct_haptic = round(haptic_pi / tm_candidates * 100, 1)

    md = textwrap.dedent(f"""\
    # CSD Scope Coverage Waterfall — CoordRep v1

    ## Summary

    CoordRep v1 converts **{valid_records:,}** CSD entries to valid records
    from the **{total:,}** entries in CSD 2024.3.  The headline figure of
    **{round(valid_records/total*100, 2)}%** of all CSD entries may give a
    misleading impression because the majority of CSD entries are purely
    organic or main-group compounds with no coordination environment.

    Using progressively narrower denominators:

    | Denominator | Entries | Valid v1 records | Coverage |
    |---|---:|---:|---:|
    | All CSD entries | {total:,} | {valid_records:,} | {round(valid_records/total*100,2)}% |
    | Transition-metal candidates | {tm_candidates:,} | {valid_records:,} | {pct_tm}% |
    | Mononuclear η1 coord. candidates | {coordination_candidates:,} | {valid_records:,} | {pct_coord}% |
    | Intended v1 domain | {intended_v1:,} | {valid_records:,} | {pct_v1}% |

    ## Hierarchical Exclusion Waterfall

    Starting from all {total:,} CSD 2024.3 entries, entries are removed
    by their **first** (highest-priority) exclusion reason:

    | Stage | Removed | Remaining |
    |---|---:|---:|
    | No 3D / atom-resolved structure | {no_3d:,} | {total - no_3d:,} |
    | No transition metal | {no_tm:,} | {tm_candidates:,} |
    | **→ TM candidates (denominator 2)** | — | **{tm_candidates:,}** |
    | Multinuclear / extended | {multinuclear:,} | {tm_candidates - multinuclear:,} |
    | Haptic / π (η > 1) | {haptic_pi:,} | {tm_candidates - multinuclear - haptic_pi:,} |
    | Disorder / partial occupancy | {disorder:,} | {coordination_candidates - disorder:,} |
    | Incompatible bonding graph | {incompatible_bonding:,} | — |
    | Other | {outside_v1_other:,} | — |
    | **→ Intended v1 domain** | — | **{intended_v1:,}** |
    | In-domain validation failure | {in_domain_failures:,} | {valid_records:,} |
    | **→ Valid CoordRep v1 records** | — | **{valid_records:,}** |

    ## Why the 8.83% is not the coverage fraction

    1. **{no_3d:,}** entries ({round(no_3d/total*100,1)}%) lack 3D coordinates entirely.
    2. **{no_tm:,}** entries ({round(no_tm/total*100,1)}%) are organic or main-group
       compounds with no transition metal.
    3. Together, these account for **{round((no_3d+no_tm)/total*100,1)}%** of the CSD —
       they were never candidates for coordination representation.

    Among the **{tm_candidates:,}** TM candidates:
    - **{pct_multi}%** are multinuclear ({multinuclear:,} entries) — outside v1 scope
      but addressed by the CoordRep-Multi-v0 prototype extension.
    - **{pct_haptic}%** have haptic/π coordination ({haptic_pi:,} entries) — outside v1
      scope but addressed by the CoordRep-Haptic-v0 prototype extension.
    - **{round(disorder/tm_candidates*100,1)}%** have crystallographic disorder
      ({disorder:,} entries).

    ## Intended v1 conversion rate

    Among the **{intended_v1:,}** entries within CoordRep v1's intended scope
    (mononuclear, η1, CN 2–6, no disorder, no polymer), **{pct_v1}%** are
    successfully converted.  The remaining **{in_domain_failures:,}**
    ({round(in_domain_failures/intended_v1*100,2)}%) fail final serialization
    validation.

    ## Files

    - `csd_scope_waterfall.csv` — Hierarchical waterfall with mutually exclusive categories
    - `csd_scope_waterfall_summary.json` — Machine-readable summary with all metrics
    - `csd_scope_exclusion_by_denominator.csv` — Each category as % of four denominators
    - `fig5B_updated_scope_waterfall.csv` — Data for revised Figure 5B
    - `csd_scope_waterfall_for_si.md` — This file
    - `README.md` — Scope notes
    """)

    (OUTDIR / "csd_scope_waterfall_for_si.md").write_text(md)
    print("  csd_scope_waterfall_for_si.md")


def write_readme():
    readme = textwrap.dedent("""\
    # CSD Scope Coverage Waterfall

    ## Purpose

    This directory provides a **multi-denominator scope coverage analysis**
    of CoordRep v1 against CSD 2024.3, responding to the reviewer question:

    > "What fraction of real coordination complexes are excluded by current
    > filtering criteria?"

    ## Key points

    - **8.83%** is the fraction of *all* CSD entries converted to valid
      CoordRep v1 records.  This is **not** the fraction of real
      coordination complexes supported.
    - The vast majority of excluded entries are organic / main-group (no TM)
      or lack 3D coordinates — they were never coordination candidates.
    - Reviewer-facing metrics should use the **TM-candidate** or
      **intended-v1-domain** denominators.

    ## Scope statement

    - CoordRep v1 is intentionally restricted to **mononuclear,
      atom-resolved η1 coordination snapshots**.
    - Multinuclear and haptic / π systems are outside v1 scope but are
      addressed by supplementary prototype extensions
      (see `coordrep_extension_prototypes/`).

    ## Exclusion priority

    When an entry has multiple exclusion reasons, it is classified by the
    **first** (highest-priority) reason encountered in the sequential
    filter pipeline:

    1. No 3D / atom-resolved structure
    2. No transition metal
    3. Multinuclear / extended coordination
    4. Haptic / π coordination (η > 1)
    5. CN outside 2–6
    6. Disorder / partial occupancy
    7. Polymeric
    8. Donor / SMILES failure
    9. Other
    10. In-domain CoordRep validation failure

    This ensures every entry appears in exactly one exclusion category.

    ## No raw CSD coordinates

    No raw CSD coordinates are exported.  All files contain aggregate
    counts and derived statistics only.

    ## Regeneration

    ```bash
    cd libcoordrep
    python scripts/generate_csd_scope_waterfall.py
    ```
    """)

    (OUTDIR / "README.md").write_text(readme)
    print("  README.md")


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("CSD Scope Coverage Waterfall")
    print(f"  Total CSD entries:  {total:,}")
    print(f"  Valid v1 records:   {valid_records:,}")
    print(f"  Whole-CSD fraction: {round(valid_records/total*100, 2)}%")
    print()

    print("Exclusion breakdown:")
    print(f"  No 3D:             {no_3d:>10,}")
    print(f"  No TM:             {no_tm:>10,}")
    print(f"  Multinuclear:      {multinuclear:>10,}")
    print(f"  Haptic/π:          {haptic_pi:>10,}")
    print(f"  Disorder:          {disorder:>10,}")
    print(f"  Incompatible:      {incompatible_bonding:>10,}")
    print(f"  Other:             {outside_v1_other:>10,}")
    print(f"  In-domain fail:    {in_domain_failures:>10,}")
    print(f"  Valid records:     {valid_records:>10,}")
    print(f"  Sum:               {no_3d+no_tm+multinuclear+haptic_pi+disorder+incompatible_bonding+outside_v1_other+in_domain_failures+valid_records:>10,}")
    print()

    print("Denominator-corrected coverage:")
    print(f"  Among TM candidates ({tm_candidates:,}):             "
          f"{round(valid_records/tm_candidates*100, 1)}%")
    print(f"  Among coord. candidates ({coordination_candidates:,}):        "
          f"{round(valid_records/coordination_candidates*100, 1)}%")
    print(f"  Among intended v1 ({intended_v1:,}):            "
          f"{round(valid_records/intended_v1*100, 1)}%")
    print()

    write_waterfall_csv()
    write_summary_json()
    write_exclusion_by_denominator()
    write_fig5b_data()
    write_si_markdown()
    write_readme()

    print(f"\nAll outputs → {OUTDIR}")
    print("Done: CSD scope waterfall")
