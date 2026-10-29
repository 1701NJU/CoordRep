# Revision Experiments

Summary of all experiments conducted during JACS revision.

## 1. CoordRep-ID L0–L3 Identity Hierarchy

| Item | Detail |
|------|--------|
| Purpose | Define multi-resolution identity keys for coordination complexes |
| Script | `scripts/benchmark_identity_layers.py`, `scripts/csd_identity_family_benchmark.py` |
| Input | tmQM + CSD retained entries |
| Output | `revision_results/csd_external/csd_identity_layer_retention.csv` |
| Key result | CSD families: L0 = 21.5%, L1 = 61.1%, L2 = 79.6%, L3 = 88.9% match |

## 2. Perturbation Robustness

| Item | Detail |
|------|--------|
| Purpose | Quantify identity key stability under coordinate noise |
| Script | `scripts/test_geometric_robustness.py` |
| Input | 200 molecules × σ = 0.001–0.2 Å |
| Output | `outputs/robustness/perturbation_sweep.csv` |
| Key result | σ = 0.01 Å: L0 = 33%, L3 = 94% retention |

## 3. Factorized Tokenizer Ablation

| Item | Detail |
|------|--------|
| Purpose | Fix composite metal token mismatch; quantify factorized improvement |
| Script | `scripts/factorized_token_ablation.py`, `scripts/tokenizer_audit.py` |
| Input | tmQM training data |
| Output | `revision_results/factorized_token/factorized_ablation_summary.csv` |
| Key result | Metal Top-1: 13.8% → 67.8%; CN: 13.8% → 49.0% |

## 4. Tool A Ablation (Field Attribution)

| Item | Detail |
|------|--------|
| Purpose | Identify which CoordRep fields drive donor prediction |
| Script | `scripts/run_tool_a_ablation.py` |
| Input | Pretrained MLM + 19,992 test positions |
| Output | `revision_results/tool_a_ablation/tool_a_ablation_summary.csv` |
| Key result | full = 85.7%, LigandFreq = 86.2%, no_lig = 61.4%, no_geom ≈ full |

## 5. Multidentate Coverage

| Item | Detail |
|------|--------|
| Purpose | Demonstrate CoordRep handles multidentate ligands |
| Script | `scripts/multidentate_coverage.py`, `scripts/multidentate_tool_a_eval.py` |
| Input | 68,288 complexes (tmQM + CSD) |
| Output | `revision_results/multidentate/multidentate_coverage_summary.csv` |
| Key result | 76.1% contain multidentate; bidentate Top-1 = 86.7% |

## 6. GNN Baselines — Donor Annotation

| Item | Detail |
|------|--------|
| Purpose | Compare CoordRep-MLM vs 2D/3D GNN for donor prediction |
| Script | `scripts/gnn_baselines/run_donor_annotation.py` |
| Input | 2000 test complexes |
| Output | `revision_results/gnn_baselines/donor_annotation_summary.csv` |
| Key result | EGNN F1 = 0.998; GIN+context = 0.904; LigandFreq = 0.896 |

## 7. GNN Baselines — Hard-Negative Ranking

| Item | Detail |
|------|--------|
| Purpose | Compare semantic ranking: GNN vs CoordRep-Ranker |
| Script | `scripts/gnn_baselines/run_hard_negative_gnn.py` |
| Input | 1501 test groups × 20 decoys |
| Output | `revision_results/gnn_baselines/3d_hard_negative_summary.csv` |
| Key result | EGNN AUROC = 0.499; CoordRep-Ranker strict = 0.800; stereo_hard = 0.948 |

## 8. CSD External Audit

| Item | Detail |
|------|--------|
| Purpose | Validate CoordRep on 200k CSD entries; quantify scope |
| Script | `scripts/csd_external_audit.py` |
| Input | 200,000 CSD entries (requires CSD license) |
| Output | `revision_results/csd_external_summary_only/summary.json` |
| Key result | 17,038 retained (8.5%); boundary fraction 13.4% |

## 9. CSD Identity Family Benchmark

| Item | Detail |
|------|--------|
| Purpose | Validate L0–L3 hierarchy on real CSD refcode families |
| Script | `scripts/csd_identity_family_benchmark.py` |
| Input | 765 CSD families with ≥2 entries |
| Output | `revision_results/csd_external_summary_only/csd_identity_layer_retention.csv` |
| Key result | L3 within-family match = 88.9% |

## 10. Tool B Grammar Repair Baselines

| Item | Detail |
|------|--------|
| Purpose | Compare rule repair vs MLM repair on same-500 subset |
| Script | `scripts/run_toolb_baselines.py`, `scripts/csd_toolb_transfer_eval.py` |
| Input | 500 synthetic + 500 CSD-derived corrupted strings |
| Output | `revision_results/toolb_sequence_baselines/toolb_same500_subset.csv` |
| Key result | Rule: valid = 1.0, field = 0.898; MLM: valid = 0.852, field = 0.606 |

## 11. CoordRep Field Ablation for Ranker

| Item | Detail |
|------|--------|
| Purpose | Prove each CoordRep field drives its corresponding decoy discrimination |
| Script | `scripts/gnn_baselines/run_hard_negative_gnn.py` (field ablation mode) |
| Input | CoordRep-Ranker + field-masked inputs |
| Output | `revision_results/gnn_baselines/taskB_coordrep_field_ablation.csv` |
| Key result | no_stereo: 0.948 → 0.500; no_metal: 0.792 → 0.500 |

## 12. CSD Casebook

| Item | Detail |
|------|--------|
| Purpose | Select and QC 4 main-text application cases |
| Script | `scripts/build_csd_casebook.py` |
| Input | CSD retained entries |
| Output | `revision_results/csd_casebook/casebook_qc_report.md` |
| Key result | ACUWOK/AFOSIA/AGOTIA/CIJWUO all QC-verified |
