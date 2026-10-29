# Manuscript Revision Map

## Abstract

| Change | Reviewer(s) |
|--------|-------------|
| Remove "InChI-like identifier" → "multi-resolution identity hierarchy" | R4.1, R2.5 |
| Remove "inverse ligand design" → "donor annotation and field attribution" | R1.3, R4.5, R2.12 |
| Add "auditable coordination-record layer" | R2.1, R2.10 |
| Add "mononuclear η1 scope" qualifier | R1.1, R2.4 |
| Add "GNN-complementary" positioning | R2.1 |

## Introduction

| Change | Reviewer(s) |
|--------|-------------|
| Clarify gap: coordination records, not merely predictive model | R2.10, R2.12 |
| Add GNN complementarity paragraph | R2.1, R2.6 |
| Add literature: Pidko, Corminboeuf, Kulik, hemilability datasets | R3.1 |
| Add literature: Fey, Cundari, geometry indices | R3.2 |
| Position relative to WL hash and graph fingerprints | R3.3, R2.11 |

## Results — Section 1: Canonicalization and Identity Hierarchy

| Change | Reviewer(s) |
|--------|-------------|
| Rename: "Canonicalization defines a multi-resolution identity hierarchy" | R4.1 |
| Insert CoordRep-ID L0–L3 definition and construction | R4.1, R2.5 |
| Add perturbation robustness: σ = 0.001–0.2 Å sweep | R2.5, R4.1 |
| Add M–L bond scaling: ±1–5% | R2.5 |
| Add CSD family benchmark: L0 = 21.5%, L3 = 88.9% | R2.5, R4.1 |
| Add ACUWOK case: 12 entries, 12 unique L0, single L3 | R2.10 |
| **Figure 2 panels c–f** | – |

## Results — Geometry and Boundary

| Change | Reviewer(s) |
|--------|-------------|
| Add boundary-aware interpretation using CShM competition | R2.3 |
| Add AFOSIA case: Δ = 0.01, boundary threshold < 1.0 | R2.3, R2.10 |
| Note: boundary cases = 13.4% of CSD retained | R2.3 |

## Results — ML / Applications

| Change | Reviewer(s) |
|--------|-------------|
| Replace Tool A "design" narrative with donor annotation + attribution | R1.3, R4.3, R4.5, R2.12 |
| Add Tool A ablation table: full/LigandFreq/no_lig/no_geom/shuffled | R1.3, R4.3 |
| Add GNN baseline comparison subsection | R2.1, R2.2 |
| — Donor annotation: EGNN F1 = 0.998 vs CoordRep ≈ 0.857 | R2.1 |
| — Hard-negative: EGNN AUROC = 0.499 vs CoordRep = 0.800 | R2.1, R2.7 |
| — CSD identity: WL/ECFP/SMILES vs CoordRep L0–L3 | R2.11, R3.3 |
| Add strict semantic decoy analysis (stereo_hard = 0.948) | R2.9 |
| Add field ablation: stereo 0.948→0.500, metal 0.792→0.500 | R2.9 |
| Add CSD application workflow: audit → casebook 4 cases | R2.10 |
| Add multidentate paragraph: 76.1% multidentate, bidentate Top-1 = 86.7% | R1.4, R4.4 |
| Add Tool B grammar repair case: AGOTIA 5/6 fields | R2.10 |
| Add stereo semantic case: CIJWUO 5/5 ranker correct | R2.10 |
| **Figure 5 panels a–h** | – |

## Results — Tokenizer

| Change | Reviewer(s) |
|--------|-------------|
| Add tokenizer audit paragraph: composite mismatch found | R1.2, R2.8 |
| Add factorized ablation: metal 67.8%, CN 49.0% | R1.2, R2.8 |
| Caveat: matched short training, not full pretrained | R1.2 |

## Methods

| Change | Reviewer(s) |
|--------|-------------|
| Tokenizer: factorized default description | R1.2 |
| CoordRep-ID construction: L0–L3 definition | R4.1 |
| GNN baselines: GIN, EGNN architecture, training details | R2.1, R2.2 |
| CSD audit: 200k pipeline, rejection criteria | R1.1, R2.4 |
| Multidentate validator: error types and rates | R1.4 |
| 3 annotated CoordRep examples (or move to Results box) | R4.2 |
| Data availability: repository URL + CSD license caveat | R3.5 |

## Discussion

| Change | Reviewer(s) |
|--------|-------------|
| v1/v2 scope table: mononuclear η1 CN 2–6 | R1.1 |
| Limitations: multinuclear, ηn, disorder, thermodynamic stability | R1.1, R2.13 |
| Spin-state caveat for 3d metals | R3.4 |
| GNN complementarity: geometry-local vs field-level semantic | R2.1, R2.6, R2.7 |
| Statistical vs physical: each field drives corresponding discrimination | R2.9 |
| Future extensions: CoordRep-Multi, Periodic, Haptic, DFT integration | R1.1, R2.13 |
| Explicitly not claiming stability prediction or electronic structure | R1.3, R2.13 |

## SI Additions

| Addition | Reviewer(s) |
|----------|-------------|
| CSD audit waterfall table with all rejection counts | R1.1, R2.4 |
| v1 scope table | R1.1 |
| Tokenizer audit detail + by-task factorized results | R1.2, R2.8 |
| Tool A ablation by CN, denticity, donor element, ligand frequency | R1.3, R4.3 |
| Multidentate coverage by CN, metal row, source; validator failures | R1.4, R4.4 |
| GNN donor annotation full table | R2.1, R2.2 |
| GNN hard-negative full table with per-decoy-type AUROC | R2.1, R2.2 |
| EGNN+metadata control experiment | R2.1 |
| Strict vs plausible decoy regrouping | R2.9 |
| CoordRep field ablation for ranker | R2.9 |
| CSD identity baselines: SMILES/WL/ECFP/CoordRep L0–L3 | R2.11, R3.3 |
| Perturbation robustness sweep table | R2.5 |
| Bond scaling table | R2.5 |
| CSD casebook: 4 main + 4 SI cases | R2.10 |
| QC report for casebook cases | R2.10 |
| 3 fully annotated CoordRep examples | R4.2 |
| Tool B repair workflow: same-500 baselines table | – |
| Ranker minphys ablation | – |
