# DATA_MANIFEST — JACS Revision Reproducibility Package

All paths are relative to `libcoordrep/`.

## Box 1: CoordRep Representative Records

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Box 1 Ex. 1–3 (v1) | `revision_results/box1_coordrep_records/box1_full_records.jsonl` | CSD-derived | selected examples | NA | NA | `scripts/export_box1_coordrep_records.py` | derived only, no coordinates | v1 production records |
| Box 1 Ex. 1–3 display | `revision_results/box1_coordrep_records/box1_display_records.md` | — | — | NA | NA | same | derived only | human-readable display |
| Box 1 Ex. 1–3 validation | `revision_results/box1_coordrep_records/box1_validation_report.csv` | — | — | NA | NA | same | derived only | per-check validation |
| Box 1 Ex. 4–5 (v2beta) | `revision_results/box1_v2beta_examples/box1_v2beta_full_records.jsonl` | CSD-derived | ETOREN, FEROCE01 | NA | NA | `scripts/generate_v2beta_postfix_audit.py` | derived only, no coordinates | v2beta extension records |
| Box 1 Ex. 4–5 display | `revision_results/box1_v2beta_examples/box1_v2beta_display_records.txt` | — | — | NA | NA | same | derived only | human-readable display |
| Box 1 Ex. 4–5 validation | `revision_results/box1_v2beta_examples/box1_v2beta_validation_report.csv` | — | — | NA | NA | same | derived only | 14-gate validation |
| Box 1 Ex. 4–5 identity | `revision_results/box1_v2beta_examples/box1_v2beta_identity_keys.csv` | — | — | NA | NA | same | derived only | L0–L3 + LocalIDs |

## Figure 2: CoordRep-ID Robustness

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 2A rigid-body | `revision_results/identity_robustness/rigid_body_invariance_test.csv` | tmQM + CSD-derived | 4000 trials | 42 | NA | `scripts/test_geometric_robustness.py` | derived only | rotation/translation invariance |
| Fig. 2B perturbation | `revision_results/identity_robustness/perturbation_sweep.csv` | tmQM + CSD-derived | perturbation sweep | 42 | NA | same | derived only | noise robustness |
| Fig. 2C bond scaling | `revision_results/identity_robustness/bond_scaling_sweep.csv` | tmQM + CSD-derived | scaling sweep | 42 | NA | same | derived only | bond-length scaling |
| Fig. 2D family benchmark | `revision_results/identity_robustness/csd_refcode_family_benchmark.csv` | CSD-derived | refcode families | NA | NA | `scripts/csd_identity_family_benchmark.py` | derived only | L0–L3 precision/recall |
| Fig. 2D family case | `revision_results/identity_robustness/family_case_study.csv` | CSD-derived | ACUWOK family | NA | NA | same | derived only | representative case |
| Fig. 2A example | `revision_results/fig2a_example/` | CSD-derived | ACUWOK | NA | NA | `scripts/_export_fig2a_example.py` | derived only | visual example |

## Figure 4: Field Learning Diagnostics

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 4B syntax validity | `revision_results/fig4_field_learning_revision/fig4B_syntax_validity_vs_steps.csv` | CoordRep train | full train | 42 | pretrain_v3 | `scripts/export_fig4_field_learning.py` | derived only | validity vs training step |
| Fig. 4B mask recovery | `revision_results/fig4_field_learning_revision/fig4B_mask_ratio_recovery.csv` | CoordRep test | test split | 42 | pretrain_v3 | same | derived only | mask ratio sweep |
| Fig. 4C donor attribution | `revision_results/fig4_field_learning_revision/fig4C_donor_field_attribution.csv` | CoordRep test | test split | 42 | pretrain_v3 | same | derived only | field ablation controls |
| Fig. 4C caption numbers | `revision_results/fig4_field_learning_revision/fig4_all_caption_numbers.json` | — | — | — | — | same | derived only | all Fig. 4 numbers |
| Fig. 4D graph baselines | `revision_results/gnn_baselines/donor_annotation_summary.csv` | CoordRep test | matched test | 2026 | coordrep_ranker | `scripts/gnn_baselines/` | derived only | E-GNN / SchNet / CoordRep comparison |
| Fig. 4D semantic decoys | `revision_results/gnn_baselines/hard_negative_summary_fixed.csv` | CoordRep test | decoy sets | 2026 | coordrep_ranker | same | derived only | strict + stereo decoys |
| Fig. 4D donor by CN | `revision_results/gnn_baselines/3d_donor_annotation_by_cn.csv` | CoordRep test | by CN | 2026 | — | same | derived only | CN-stratified results |
| Fig. 4D donor by dent. | `revision_results/gnn_baselines/3d_donor_annotation_by_denticity.csv` | CoordRep test | by denticity | 2026 | — | same | derived only | denticity-stratified |

## Table 2: Tokenizer Ablation

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Table 2 main | `revision_results/factorized_token/factorized_ablation_summary.csv` | CoordRep test | test split | 42 | pretrain_v3 + ablations | `scripts/factorized_token_ablation.py` | derived only | composite vs factorized |
| Table 2 by task | `revision_results/factorized_token/factorized_by_task.csv` | — | — | 42 | — | same | derived only | per-task breakdown |
| Table 2 by CN | `revision_results/factorized_token/factorized_by_cn.csv` | — | — | 42 | — | same | derived only | per-CN breakdown |
| Table 2 by denticity | `revision_results/factorized_token/factorized_by_denticity.csv` | — | — | 42 | — | same | derived only | per-denticity |
| Table 2 summary | `revision_results/factorized_token/summary.json` | — | — | — | — | same | derived only | aggregate numbers |
| Table 2 comparison | `revision_results/factorized_token/tool_a_standard_eval_comparison.csv` | — | — | 42 | pretrain_v3 | `scripts/factorized_token_qc_and_standard_eval.py` | derived only | old vs new tokenizer |

## Table 3: Multidentate Donor-Set Evaluation

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Table 3 main | `revision_results/multidentate/tool_a_by_denticity_final.csv` | CoordRep test | multidentate subset | 42 | pretrain_v3 | `scripts/multidentate_tool_a_eval.py` | derived only | Tool A by denticity |
| Table 3 coverage | `revision_results/multidentate/multidentate_coverage_summary.csv` | CSD-derived | full CSD | NA | NA | `scripts/multidentate_coverage.py` | derived only | coverage by denticity |
| Table 3 invalid cases | `revision_results/multidentate/tool_a_multidentate_invalid_cases.csv` | CoordRep test | error cases | 42 | pretrain_v3 | same | derived only | hard error examples |
| Table 3 validator | `revision_results/multidentate/multidentate_validator_summary.csv` | CoordRep test | — | NA | NA | `scripts/multidentate_validator.py` | derived only | validator pass rates |

## Figure 5B: Full-CSD Scope Audit

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 5B v1 scope | `revision_results/csd_scope_waterfall_revision/csd_scope_waterfall.csv` | CSD 2024.3 | 1,413,222 entries | NA | NA | `scripts/generate_csd_scope_waterfall.py` | derived only | v1 waterfall: 124,837 records |
| Fig. 5B scope summary | `revision_results/csd_scope_waterfall_revision/csd_scope_waterfall_summary.json` | — | — | NA | NA | same | derived only | v1 8.83% coverage |
| Fig. 5B merged panel | `revision_results/figure_ready_merged_fig5/fig5B_full_csd_scope_audit.csv` | — | — | NA | NA | `scripts/build_merged_fig5.py` | derived only | figure-ready CSV |
| Fig. 5B merged JSON | `revision_results/figure_ready_merged_fig5/fig5B_full_csd_scope_audit.json` | — | — | NA | NA | same | derived only | key numbers |

### Key numbers (Fig. 5B)
- v1 valid records: **124,837**
- All-CSD coverage: **8.83%**

## Figure 5C: Boundary Geometry Atlas

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 5C CN4 atlas | `revision_results/csd_pathfinder_full/cn4_sp_td_atlas.csv` | CSD-derived | CN=4 | NA | NA | `scripts/run_csd_pathfinder.py` | derived only | 47,187 records |
| Fig. 5C CN5 atlas | `revision_results/csd_pathfinder_full/cn5_tbpy_spy_atlas.csv` | CSD-derived | CN=5 | NA | NA | same | derived only | 14,897 records |
| Fig. 5C CN6 atlas | `revision_results/csd_pathfinder_full/cn6_oh_distortion_atlas.csv` | CSD-derived | CN=6 | NA | NA | same | derived only | 41,381 records |
| Fig. 5C ridge enrichment | `revision_results/csd_pathfinder_full/cn5_pathway_ridge_summary.csv` | — | CN=5 ridge | NA | NA | same | derived only | 42.2% intermediate |

### Key numbers (Fig. 5C)
- Boundary records: **17,791**
- Boundary fraction: **14.3%**
- CN5 boundary/ridge fraction: **42.2%**

## Figure 5D: L3 Family Geometry Trajectories

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 5D trajectories | `revision_results/csd_pathfinder_full/l3_family_geometry_trajectories.csv` | CSD-derived | L3 families | NA | NA | `scripts/export_fig5d_family_trajectory.py` | derived only | 6,358 families |
| Fig. 5D merged panel | `revision_results/figure_ready_merged_fig5/fig5D_family_polymorphism_summary.csv` | — | — | NA | NA | `scripts/build_merged_fig5.py` | derived only | figure-ready |

### Key numbers (Fig. 5D)
- Nontrivial L3 families: **6,358**

## Figure 5E: v2-Beta Extension Audit

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 5E postfix audit | `revision_results/full_csd_v2beta_postfix_audit/multinuclear_molecular_subset_audit.csv` | CSD-derived | 4,579 molecular multinuclear | NA | NA | `scripts/generate_v2beta_postfix_audit.py` | derived only | 100% conversion/roundtrip |
| Fig. 5E scope reclass. | `revision_results/full_csd_v2beta_postfix_audit/multinuclear_scope_reclassification.csv` | CSD-derived | 5,000 sample | NA | NA | same | derived only | 91.6% molecular |
| Fig. 5E coverage | `revision_results/full_csd_v2beta_postfix_audit/v2beta_postfix_coverage_summary.csv` | CSD-derived | full CSD | NA | NA | same | derived only | 500,299 total |
| Fig. 5E claim | `revision_results/full_csd_v2beta_postfix_audit/v2beta_postfix_claim_recommendation.json` | — | — | NA | NA | same | derived only | claim level C |
| Fig. 5E haptic audit | `revision_results/full_csd_v2beta_audit/coordrep_haptic_full_audit.csv` | CSD-derived | 2,000 sample | NA | NA | `scripts/generate_full_csd_v2beta_audit.py` | derived only | haptic validation |

### Key numbers (Fig. 5E)
- v1 + v2beta records: **500,299**
- All-CSD coverage: **8.83% → 35.40%**
- TM candidate coverage: **20.3% → 81.28%**
- Coverage gain: **+300.8%**

## Figure 5G: CoordRep-Rosetta Hard-Control Benchmark

| Manuscript item | File path | Source dataset | Subset/split | Seed | Checkpoint | Script | CSD redistribution status | Notes |
|---|---|---|---|---|---|---|---|---|
| Fig. 5G hard control | `revision_results/coordrep_rosetta_hard_controls/hard_control_summary.json` | CSD-derived | 143 hard pools | 42 | NA | `scripts/generate_rosetta_hard_controls.py` | derived only | GO decision |
| Fig. 5G method comparison | `revision_results/coordrep_rosetta_hard_controls/hard_control_method_comparison.csv` | — | — | 42 | NA | same | derived only | 5 methods |
| Fig. 5G pool index | `revision_results/coordrep_rosetta_hard_controls/hard_pool_index.csv` | CSD-derived | — | 42 | NA | same | derived only | pool definitions |
| Fig. 5G detailed results | `revision_results/coordrep_rosetta_hard_controls/hard_control_detailed_results.csv` | — | — | 42 | NA | same | derived only | per-pool results |
| Fig. 5G manual audit | `revision_results/coordrep_rosetta_hard_controls/hard_control_manual_audit.csv` | CSD-derived | 90 sampled pairs | 42 | NA | same | derived only | chemical homology |
| Fig. 5G audit summary | `revision_results/coordrep_rosetta_hard_controls/hard_control_manual_audit_summary.json` | — | — | — | NA | same | derived only | audit rates |
| Fig. 5G figure bars | `revision_results/coordrep_rosetta_hard_controls/fig_hard_control_bars.csv` | — | — | — | NA | same | derived only | figure-ready |
| Fig. 5G drawing pack | `revision_results/fig5F_drawing_pack/fig5F_benchmark_summary.csv` | — | — | — | NA | — | derived only | for figure designer |
| Fig. 5G pilot | `revision_results/coordrep_rosetta_active_site_homology/coordrep_rosetta_summary.json` | CSD-derived | 500 MOF nodes | 42 | NA | `scripts/generate_rosetta_homology.py` | derived only | SI proof-of-concept |

### Key numbers (Fig. 5G)
- Full CoordRep Top-1: **86.7%**
- Full CoordRep AUROC: **0.949**
- Top-1 gain over donor-only baseline: **+18.9 pp**
- Manual audit, CoordRep Top-1 homologous: **93.3%**
- Manual audit, donor baseline homologous: **56.7%**
- Manual audit, hard negatives homologous: **0.0%**

## Supplementary Figures

| Manuscript item | File path | Source dataset | Script | CSD redistribution status |
|---|---|---|---|---|
| Fig. S2 | `revision_results/si_figures/Fig_S2.pdf` | CSD-derived | `scripts/si_figures/` | derived only |
| Fig. S3 | `revision_results/si_figures/Fig_S3.pdf` | CSD-derived | `scripts/si_figures/` | derived only |
| Fig. S4 | `revision_results/si_figures/Fig_S4.pdf` | CSD-derived | `scripts/si_figures/` | derived only |
| Fig. S5 | `revision_results/si_figures/Fig_S5.pdf` | CSD-derived | `scripts/si_figures/` | derived only |
| Fig. S6 | `revision_results/si_figures/Fig_S6.pdf` | CSD-derived | `scripts/si_figures/` | derived only |
| Fig. S7 | `revision_results/si_figures/Fig_S7.pdf` | CSD-derived | `scripts/si_figures/` | derived only |

## Checkpoints

| Item | File path | Notes |
|---|---|---|
| Pretrained MLM | `checkpoints/pretrain_v3/config.json` | Architecture config |
| Tokenizer | `checkpoints/pretrain_v3/tokenizer.json` | Factorized tokenizer vocab |

## CSD Redistribution Status

All files in this repository contain **derived data only**:
- CoordRep string records (not coordinates)
- Hashed identity keys (irreversible)
- Aggregate statistics
- Validation labels
- CShM descriptors
- Retrieval scores

See `CSD_REDISTRIBUTION_NOTICE.md` for full details.
