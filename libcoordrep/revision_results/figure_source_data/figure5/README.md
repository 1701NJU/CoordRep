# Figure 5: Full-CSD Audit, Boundary Geometry, and Rosetta

Source data files for Figure 5 panels.

## 5B: Scope Audit

| File | Description |
|---|---|
| `../../csd_scope_waterfall_revision/csd_scope_waterfall.csv` | v1 scope waterfall |
| `../../csd_scope_waterfall_revision/csd_scope_waterfall_summary.json` | v1 summary (124,837 records, 8.83%) |
| `../../figure_ready_merged_fig5/fig5B_full_csd_scope_audit.csv` | Figure-ready CSV |

## 5C: Boundary Geometry Atlas

| File | Description |
|---|---|
| `../../csd_pathfinder_full/cn4_sp_td_atlas.csv` | CN4 SP↔Td atlas (47,187) |
| `../../csd_pathfinder_full/cn5_tbpy_spy_atlas.csv` | CN5 TBPy↔SPY atlas (14,897) |
| `../../csd_pathfinder_full/cn6_oh_distortion_atlas.csv` | CN6 Oh distortion atlas (41,381) |
| `../../csd_pathfinder_full/cn5_pathway_ridge_summary.csv` | CN5 ridge enrichment (42.2%) |

## 5D: L3 Family Trajectories

| File | Description |
|---|---|
| `../../csd_pathfinder_full/l3_family_geometry_trajectories.csv` | 6,358 L3 families |
| `../../figure_ready_merged_fig5/fig5D_family_polymorphism_summary.csv` | Figure-ready summary |

## 5E: v2-Beta Extension

| File | Description |
|---|---|
| `../../full_csd_v2beta_postfix_audit/v2beta_postfix_coverage_summary.csv` | 500,299 total; 35.40% coverage |
| `../../full_csd_v2beta_postfix_audit/multinuclear_molecular_subset_audit.csv` | 4,579 molecular multinuclear |
| `../../full_csd_v2beta_postfix_audit/multinuclear_scope_reclassification.csv` | Scope reclassification |

## 5F: Boundary Competition

| File | Description |
|---|---|
| `../../figure_ready_merged_fig5/fig5F_boundary_competition_summary.json` | Boundary statistics |
| `../../figure_ready_merged_fig5/fig5F_quantitative_summary_table.csv` | Quantitative summary |

## 5G: CoordRep-Rosetta

| File | Description |
|---|---|
| `../../coordrep_rosetta_hard_controls/hard_control_summary.json` | Hard-control benchmark (GO) |
| `../../coordrep_rosetta_hard_controls/hard_control_method_comparison.csv` | 5 methods comparison |
| `../../coordrep_rosetta_hard_controls/hard_control_manual_audit.csv` | Manual audit (90 pairs) |
| `../../coordrep_rosetta_hard_controls/fig_hard_control_bars.csv` | Figure-ready bars |
| `../../fig5F_drawing_pack/fig5F_benchmark_summary.csv` | Drawing pack for designer |

### Key numbers
- Boundary records: **17,791** (14.3%)
- CN5 boundary/ridge: **42.2%**
- Nontrivial L3 families: **6,358**
- v1+v2beta total: **500,299** (35.40% all-CSD)
- Rosetta top-1: **86.7%**, AUROC: **0.949**, delta: **+18.9 pp**
