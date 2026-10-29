#!/usr/bin/env python3
"""
CoordRep JACS Revision Package — Integrity Checker

Verifies that all required source-data files exist, no raw CSD files are
present, and manifest paths are valid.

Usage:
    python scripts/check_revision_package.py
"""

import json
import os
import sys
from pathlib import Path

# Resolve base directory (libcoordrep/)
BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "revision_results"

PASS = 0
FAIL = 0


def check(label: str, condition: bool, detail: str = ""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {label}")
    else:
        FAIL += 1
        msg = f"  FAIL: {label}"
        if detail:
            msg += f" — {detail}"
        print(msg)


def check_file(label: str, path: Path):
    try:
        rel = path.relative_to(BASE)
    except ValueError:
        rel = path.name
    check(label, path.exists() and path.stat().st_size > 0,
          f"missing or empty: {rel}")


# ════════════════════════════════════════════════════════════════════
print("=" * 60)
print("CoordRep JACS Revision Package — Integrity Check")
print("=" * 60)

# ── 1. Box 1 records ──
print("\n[1] Box 1 records")
check_file("Box 1 v1 records (JSONL)",
           RESULTS / "box1_coordrep_records" / "box1_full_records.jsonl")
check_file("Box 1 v1 validation",
           RESULTS / "box1_coordrep_records" / "box1_validation_report.csv")
check_file("Box 1 v1 display",
           RESULTS / "box1_coordrep_records" / "box1_display_records.md")
check_file("Box 1 v2beta records (JSONL)",
           RESULTS / "box1_v2beta_examples" / "box1_v2beta_full_records.jsonl")
check_file("Box 1 v2beta validation",
           RESULTS / "box1_v2beta_examples" / "box1_v2beta_validation_report.csv")
check_file("Box 1 v2beta identity keys",
           RESULTS / "box1_v2beta_examples" / "box1_v2beta_identity_keys.csv")
check_file("Box 1 v2beta display",
           RESULTS / "box1_v2beta_examples" / "box1_v2beta_display_records.txt")

# ── 2. Figure 2 ──
print("\n[2] Figure 2 source data")
fig2_dir = RESULTS / "identity_robustness"
check_file("Rigid body invariance",
           fig2_dir / "rigid_body_invariance_test.csv")
check_file("Perturbation sweep",
           fig2_dir / "perturbation_sweep.csv")
check_file("Bond scaling sweep",
           fig2_dir / "bond_scaling_sweep.csv")
check_file("Refcode family benchmark",
           fig2_dir / "csd_refcode_family_benchmark.csv")

# ── 3. Figure 4 ──
print("\n[3] Figure 4 source data")
fig4_dir = RESULTS / "fig4_field_learning_revision"
check_file("Fig 4B syntax validity",
           fig4_dir / "fig4B_syntax_validity_vs_steps.csv")
check_file("Fig 4B mask recovery",
           fig4_dir / "fig4B_mask_ratio_recovery.csv")
check_file("Fig 4C donor attribution",
           fig4_dir / "fig4C_donor_field_attribution.csv")
check_file("Fig 4 caption numbers",
           fig4_dir / "fig4_all_caption_numbers.json")

gnn_dir = RESULTS / "gnn_baselines"
check_file("Fig 4D donor annotation summary",
           gnn_dir / "donor_annotation_summary.csv")
check_file("Fig 4D semantic decoys",
           gnn_dir / "hard_negative_summary_fixed.csv")

# ── 4. Table 2 ──
print("\n[4] Table 2 source data")
tok_dir = RESULTS / "factorized_token"
check_file("Factorized ablation summary",
           tok_dir / "factorized_ablation_summary.csv")
check_file("Factorized by task",
           tok_dir / "factorized_by_task.csv")
check_file("Token summary JSON",
           tok_dir / "summary.json")

# ── 5. Table 3 ──
print("\n[5] Table 3 source data")
multi_dir = RESULTS / "multidentate"
check_file("Tool A by denticity",
           multi_dir / "tool_a_by_denticity_final.csv")
check_file("Multidentate coverage",
           multi_dir / "multidentate_coverage_summary.csv")
check_file("Invalid cases",
           multi_dir / "tool_a_multidentate_invalid_cases.csv")

# ── 6. Figure 5B ──
print("\n[6] Figure 5B scope audit")
scope_dir = RESULTS / "csd_scope_waterfall_revision"
check_file("Scope waterfall CSV",
           scope_dir / "csd_scope_waterfall.csv")
check_file("Scope waterfall summary",
           scope_dir / "csd_scope_waterfall_summary.json")

# ── 7. Figure 5C–D ──
print("\n[7] Figure 5C–D boundary/trajectory")
path_dir = RESULTS / "csd_pathfinder_full"
check_file("CN4 SP/Td atlas",
           path_dir / "cn4_sp_td_atlas.csv")
check_file("CN5 TBPy/SPY atlas",
           path_dir / "cn5_tbpy_spy_atlas.csv")
check_file("CN6 Oh atlas",
           path_dir / "cn6_oh_distortion_atlas.csv")
check_file("L3 family trajectories",
           path_dir / "l3_family_geometry_trajectories.csv")
check_file("CN5 ridge summary",
           path_dir / "cn5_pathway_ridge_summary.csv")

# ── 8. Figure 5E ──
print("\n[8] Figure 5E v2beta audit")
v2_dir = RESULTS / "full_csd_v2beta_postfix_audit"
check_file("Multinuclear audit",
           v2_dir / "multinuclear_molecular_subset_audit.csv")
check_file("Scope reclassification",
           v2_dir / "multinuclear_scope_reclassification.csv")
check_file("Coverage summary",
           v2_dir / "v2beta_postfix_coverage_summary.csv")
check_file("Claim recommendation",
           v2_dir / "v2beta_postfix_claim_recommendation.json")

# ── 9. Figure 5G Rosetta ──
print("\n[9] Figure 5G Rosetta")
ros_dir = RESULTS / "coordrep_rosetta_hard_controls"
check_file("Hard control summary",
           ros_dir / "hard_control_summary.json")
check_file("Method comparison",
           ros_dir / "hard_control_method_comparison.csv")
check_file("Pool index",
           ros_dir / "hard_pool_index.csv")
check_file("Detailed results",
           ros_dir / "hard_control_detailed_results.csv")
check_file("Manual audit",
           ros_dir / "hard_control_manual_audit.csv")
check_file("Audit summary",
           ros_dir / "hard_control_manual_audit_summary.json")

# ── 10. Supplementary figures ──
print("\n[10] Supplementary figures")
si_dir = RESULTS / "si_figures"
for i in range(2, 8):
    check_file(f"Fig S{i}",
               si_dir / f"Fig_S{i}.pdf")

# ── 11. Merged figure 5 ──
print("\n[11] Merged Figure 5 panels")
merged = RESULTS / "figure_ready_merged_fig5"
check_file("Fig5B merged",
           merged / "fig5B_full_csd_scope_audit.csv")
check_file("Fig5D merged",
           merged / "fig5D_family_polymorphism_summary.csv")
check_file("Fig5 caption numbers",
           merged / "fig5_all_caption_numbers.json")

# ── 12. No raw CSD files ──
print("\n[12] CSD redistribution compliance")
raw_exts = {".cif", ".hkl", ".res", ".fcf"}
raw_files = []
for root, dirs, files in os.walk(BASE):
    # Skip .git
    dirs[:] = [d for d in dirs if d != ".git"]
    for f in files:
        if Path(f).suffix.lower() in raw_exts:
            raw_files.append(os.path.join(root, f))
check("No raw CSD files (*.cif, *.hkl, *.res, *.fcf)",
      len(raw_files) == 0,
      f"found: {raw_files[:5]}")

# ── 13. Key documentation ──
print("\n[13] Documentation")
repo_root = BASE.parent
check_file("README.md", repo_root / "README.md")
check_file("DATA_MANIFEST.md", repo_root / "DATA_MANIFEST.md")
check_file("CSD_REDISTRIBUTION_NOTICE.md",
           repo_root / "CSD_REDISTRIBUTION_NOTICE.md")
check_file("CODE_MAP.md", repo_root / "CODE_MAP.md")

# ── 14. Key manuscript numbers spot-check ──
print("\n[14] Key number verification")
try:
    with open(merged / "fig5_all_caption_numbers.json") as f:
        nums = json.load(f)
    check("v1 records = 124,837",
          nums.get("valid_coordrep_records") == 124837)
    check("Boundary records = 17,791",
          nums.get("boundary_records") == 17791)
    check("Boundary fraction = 14.3%",
          abs(nums.get("boundary_fraction_percent", 0) - 14.3) < 0.1)
    check("CN5 ridge = 42.2%",
          abs(nums.get("cn5_intermediate_fraction_percent", 0) - 42.2) < 0.1)
    check("L3 families = 6,358",
          nums.get("nontrivial_L3_families") == 6358)
except Exception as e:
    check("Fig5 numbers JSON readable", False, str(e))

try:
    with open(ros_dir / "hard_control_summary.json") as f:
        ros = json.load(f)
    ros_top1 = ros.get("methods_best_subtask", {}).get(
        "coordrep_full", {}).get("top1_correct_rate", 0)
    check("Rosetta top-1 = 86.7%",
          abs(ros_top1 * 100 - 86.7) < 0.5)
except Exception as e:
    check("Rosetta numbers JSON readable", False, str(e))

# ════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print(f"Results: {PASS} passed, {FAIL} failed")
print("=" * 60)

if FAIL > 0:
    print("\n⚠  Some checks failed. Review output above.")
    sys.exit(1)
else:
    print("\n✓  All checks passed. Package is ready for submission.")
    sys.exit(0)
