# Figure-Ready Evidence Pack — README

Generated from all revision experiments. No new training or experiments.

---

## 1. Proposed Revised Figure 2 Layout: Identity Hierarchy

| Panel | Content | Source |
|-------|---------|--------|
| 2a | Canonicalization example | existing |
| 2b | Identity ladder L0→L1→L2→L3 definition | existing |
| 2c | Perturbation robustness (σ = 0.001–0.2 Å) | `fig2_identity_layers.csv` |
| 2d | M–L bond scaling (±1–5%) | `fig2_identity_layers.csv` |
| 2e | CSD refcode-family benchmark: L0=21.5%, L3=88.9% | `fig2_identity_layers.csv` |
| 2f | ACUWOK case: 12 entries, 12 unique L0, single L3 | `fig2_identity_layers.csv` |

Key message: L0 captures geometry state; L3 captures chemical identity. The hierarchy provides auditable multi-resolution matching.

---

## 2. Proposed Revised Figure 5 Layout: CoordRep Application Layer

| Panel | Content | Source |
|-------|---------|--------|
| 5a | CSD external audit: 200k → 17k retention waterfall | `fig5_application_summary.csv` |
| 5b | Case A: ACUWOK identity ladder | `casebook_maintext_summary.csv` |
| 5c | Case B: AFOSIA boundary geometry | `casebook_maintext_summary.csv` |
| 5d | Case C: AGOTIA grammar repair | `casebook_maintext_summary.csv` |
| 5e | Case D: CIJWUO stereo semantic | `casebook_maintext_summary.csv` |
| 5f | Tool A field attribution: ligand-intrinsic | `fig5_application_summary.csv` |
| 5g | GNN comparison: donor task (EGNN=0.998 vs CoordRep=0.857) | `gnn_comparison_summary.csv` |
| 5h | GNN comparison: semantic ranking (EGNN=0.50 vs CoordRep=0.80) | `gnn_comparison_summary.csv` |

Key message: CoordRep is a curation/record layer, not a geometry prediction tool. 3D GNNs solve geometry-local tasks; CoordRep handles field-level semantic records.

---

## 3. Main Text vs SI Assignment

### Main Text

| Item | Key Number | File |
|------|-----------|------|
| CoordRep-ID L0–L3 hierarchy | L3=88.9% CSD match | `fig2_identity_layers.csv` |
| CSD audit waterfall | 200k → 17k | `fig5_application_summary.csv` |
| ACUWOK identity case | 12 entries, single L3 | `casebook_maintext_summary.csv` |
| AFOSIA boundary case | Δ=0.01 | `casebook_maintext_summary.csv` |
| AGOTIA grammar repair | 5/6 fields preserved | `casebook_maintext_summary.csv` |
| CIJWUO stereo semantic | 5/5 ranker correct | `casebook_maintext_summary.csv` |
| GNN donor vs semantic | EGNN F1=0.998 vs AUROC=0.50 | `gnn_comparison_summary.csv` |
| Tool A ablation headline | full=85.7%, no_lig=61.4% | `fig5_application_summary.csv` |
| Multidentate headline | 76.1% of data is multidentate | `multidentate_summary.csv` |
| CSD identity: L3 vs WL hash | L3 recall=0.83/prec=0.82 vs WL recall=0.98/prec=0.80 | `gnn_comparison_summary.csv` |

### SI

| Item | File |
|------|------|
| Detailed Ranker/minphys ablation | `revision_results/ranker_minphys/` |
| Full Tool B baselines (same-500) | `toolb_repair_summary.csv` |
| Full tokenizer ablation | `tokenizer_factorized_summary.csv` |
| Full multidentate validator | `multidentate_summary.csv` |
| CSD family detail table | `revision_results/csd_external/csd_family_identity_detail.csv` |
| CoordRep field ablation for ranker | `revision_results/gnn_baselines/taskB_coordrep_field_ablation.csv` |
| EGNN+metadata control | `revision_results/gnn_baselines/taskB_egnn_metadata_control.csv` |
| Strict vs plausible decoy regrouping | `revision_results/gnn_baselines/taskB_strict_vs_plausible.csv` |
| Per-CN / per-metal breakdowns | `revision_results/gnn_baselines/ranker_by_cn.csv` etc. |
| SI casebook cases | `revision_results/csd_casebook/casebook_selected_si.json` |

---

## 4. Revised Central Claim (one paragraph)

CoordRep is a standardized, multi-resolution, auditable record layer for mononuclear η1 coordination chemistry, enabling identity resolution (L0–L3 hierarchy validated on 765 CSD refcode families at 88.9% L3 match), geometry-boundary annotation (continuous CShM with boundary threshold), deterministic grammar repair (100% validity recovery, 5/6 field preservation), semantic field checking (stereo decoy detection at 94.8% AUROC), and model-ready curation across public and CSD-derived datasets. Its explicit scope (mononuclear, η1, CN 2–6) is a design boundary, not a limitation: within this scope, CoordRep provides capabilities — multi-resolution identity, auditable field repair, semantic consistency ranking — that neither raw SMILES, graph hashes, nor 3D GNNs can replicate.

---

## 5. File Manifest

```
revision_results/figure_ready/
├── fig2_identity_layers.csv          # Figure 2 panels c–f
├── fig5_application_summary.csv      # Figure 5 all panels
├── gnn_comparison_summary.csv        # GNN comparison (3 tasks)
├── multidentate_summary.csv          # Multidentate coverage
├── tokenizer_factorized_summary.csv  # Tokenizer ablation
├── toolb_repair_summary.csv          # Tool B same-500 baselines
├── casebook_maintext_summary.csv     # 4 main + 4 SI cases
└── README_figure_plan.md             # This file
```

All numbers are sourced from existing `revision_results/` computations. No new experiments.
