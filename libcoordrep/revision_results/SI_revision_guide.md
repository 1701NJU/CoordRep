# Supplementary Information — Revision Structure & Content Guide

This document maps the new SI outline to available data, specifies what content goes in each section, and identifies which figures/tables to include.

---

## Existing SI (pre-revision) → New SI mapping

| Old SI section | New SI section | Notes |
|---|---|---|
| SM1 (Data provenance) | Supp. Methods 1 | Expand with CSD scan waterfall |
| SM2 (Conversion pipeline) | Supp. Methods 2 | Add Box 1 records, full serializer outputs |
| SM3 (CShM) | Supp. Methods 5 | Renumber; add PathFinder details |
| SM4 (Canonicalization) | Supp. Methods 3 | Move up; add CoordRep-ID hash construction |
| SM5 (Grammar/tokenizer) | Supp. Methods 2 | Merge into Methods 2 |
| SM6 (Pretraining tasks) | Supp. Methods 6 + 7 | Split into masked-field + donor-field |
| SR1 (Corruption suite) | Supp. Results 8 | Rename: syntax repair casebook |
| SR2 (Baselines/ablations) | Supp. Results 5 + 7 | Split by topic |
| SR3 (Training setup) | Supp. Methods 6 (training) | Move to methods |
| SR4 (Reproducibility) | Supp. Methods 11 | Renumber |

---

## New Section-by-Section Content Plan

---

### Supplementary Methods 1. Data sources, scope definition, and filtering rules

**Purpose**: Document every filtering step from raw databases to the final CoordRep corpus.

**Content to include**:
1. Source databases and version dates (CSD 2024, COD snapshot, tmQM)
2. Scope definition table → `revision_results/scope/coordrep_v1_scope_table.csv`
3. Waterfall filter: full-CSD scan → retained entries
   - Data: `revision_results/csd_pathfinder_full/full_csd_filtering_waterfall.csv`
   - Summary: `revision_results/csd_pathfinder_full/full_csd_scan_summary.json`
4. Rejection reasons breakdown → `revision_results/csd_external_summary_only/csd_rejection_reasons.csv`
5. Stratified split statistics (existing Tables S1, S2 — keep)

**Tables**:
- **Table S1**: Stratified split by metal series (keep existing)
- **Table S2**: Stratified split by CN (keep existing)
- **NEW Table S3**: Full-CSD waterfall filter stages and counts (from `full_csd_filtering_waterfall.csv`)

**Figures**: None needed (data is tabular).

---

### Supplementary Methods 2. CoordRep grammar, tokenization, and complete serializer outputs

**Purpose**: Full grammar spec + Box 1 complete records (reviewer request).

**Content to include**:
1. Token types table (existing Table S4 — keep, renumber)
2. BNF grammar (keep existing SM5.2)
3. Vocabulary construction (keep existing SM5.3)
4. Composite token design rationale (keep existing SM5.4)
5. **NEW**: Box 1 full serializer outputs (SI version)
   - All 4 records from `revision_results/box1_coordrep_records/box1_si_full_records.txt`
   - Include full unwrapped strings, L0–L3 hashes
   - Reference: `box1_full_records.jsonl` for machine-readable version
6. Tokenizer implementation notes (keep existing SM5.6–5.7)

**Tables**:
- **Table S4**: Token Types (renumber from existing)

**Box/Insert**:
- SI Box: Full machine-readable CoordRep records for Box 1 A–D
  (content from `box1_si_full_records.txt`)

---

### Supplementary Methods 3. Canonicalization and CoordRep-ID construction

**Purpose**: Canonicalization algorithm + how L0/L1/L2/L3 identity keys are computed.

**Content to include**:
1. Keep existing SM4.1–4.8 (determinism guarantees, ligand ordering, etc.)
2. **NEW section SM3.x**: CoordRep-ID hierarchy construction
   - L0 = full canonical string (StateKey)
   - L1 = ShapeID (metal + shape bin + topology hash)
   - L2 = TopoID (metal + CN + shape + sorted ligand SMILES)
   - L3 = ConnID (metal + CN + sorted ligand molecular formulas)
   - Hash algorithm, truncation, and collision analysis
3. Reference → `coordrep/identity/identity_keys.py` for pseudocode
4. Examples from Box 1 showing L0–L3 at each level

**Tables**:
- **Table S5**: CoordRep-ID level definitions (new — from identity_keys module)

---

### Supplementary Methods 4. Coordination Identity Challenge

**Purpose**: Describe the representation-gap benchmark design.

**Content to include**:
1. Motivation: what operations should a "complete" coordination representation support?
2. Operation taxonomy (8 operations from `representation_capability_matrix.csv`)
3. Challenge pair construction methodology
4. Evaluation protocol: which representations can distinguish/fail each operation
5. Ground-truth source and evaluation criteria

**Data sources**:
- `revision_results/representation_gap_benchmark/representation_capability_matrix.csv`
- `revision_results/representation_gap_benchmark/coordination_identity_challenge_summary.csv`
- `revision_results/representation_gap_benchmark/challenge_pair_examples.jsonl`
- `revision_results/representation_gap_benchmark/fig_representation_gap_heatmap.csv`

**Tables**:
- **Table S6**: Full representation capability matrix (expanded version of Fig. 2 heatmap)

---

### Supplementary Methods 5. Continuous geometry and CShM analyses

**Purpose**: Full CShM methodology (keep existing SM3 content) + PathFinder details.

**Content to include**:
1. CShM mathematical definition (keep existing SM3.2–3.5)
2. Reference polyhedra library (keep Table S3 → renumber as S7)
3. **NEW**: PathFinder methodology
   - CN4 SP↔Td pathways
   - CN5 SPY↔TBPY Berry pathway
   - CN6 Oh distortion atlas
   - Boundary detection and classification
4. Delta thresholds and boundary-flag logic

**Data sources**:
- `revision_results/csd_pathfinder_full/cn4_sp_td_atlas.csv`
- `revision_results/csd_pathfinder_full/cn5_tbpy_spy_atlas.csv`
- `revision_results/csd_pathfinder_full/cn5_pathway_ridge_summary.csv`
- `revision_results/csd_pathfinder_full/cn6_oh_distortion_atlas.csv`
- `revision_results/csd_pathfinder_full/geometry_region_enrichment_cn*.csv`

**Tables**:
- **Table S7**: Reference polyhedra library (existing S3, renumbered)
- **Table S8**: PathFinder ridge-point statistics for CN5

**Figures**:
- **Fig. S1**: CN4 SP↔Td atlas (hexbin density from `cn4_sp_td_atlas.csv`)
- **Fig. S2**: CN5 Berry pathway ridge (from `cn5_pathway_ridge_points.csv`)
- **Fig. S3**: CN6 Oh distortion atlas (from `cn6_oh_distortion_atlas.csv`)

---

### Supplementary Methods 6. Masked-field learning and tokenizer ablations

**Purpose**: Full training protocol + factorized tokenizer ablation methodology.

**Content to include**:
1. Pretraining objective: MLM on CoordRep strings (existing SM6.1)
2. Task suite: three masking modes (existing SM6.2)
3. Dynamic masking and reproducibility (existing SM6.3)
4. **NEW**: Factorized vs composite tokenizer design
   - Motivation: independent recovery of metal and CN tokens
   - Training details for factorized variant
   - Data: `revision_results/factorized_token/factorized_ablation_summary.csv`
5. Training hyperparameters (existing SR3 → move here)
6. Convergence behavior (existing SR3.4)

**Tables**:
- **Table S9**: Hyperparameters (existing S8, renumbered)
- **Table S10**: Training Dynamics (existing S9, renumbered)
- **NEW Table S11**: Factorized tokenizer ablation results (from `factorized_ablation_summary.csv`)

**Figures**:
- **Fig. S4**: Pretraining dynamics (existing S1, renumbered)

---

### Supplementary Methods 7. Donor-field attribution and multidentate controls

**Purpose**: Full methodology for donor recovery evaluation and multidentate analysis.

**Content to include**:
1. Donor-field attribution protocol (masking strategies from Fig. 4C)
2. Ligand-dependence control (ligand-masked condition)
3. Multidentate donor-set recovery protocol
   - Hard-error definition (chemically invalid completions)
   - Coverage statistics: `revision_results/multidentate/multidentate_coverage_summary.csv`
4. Evaluation metrics: top-1, top-5, hard-error rate
5. Baseline definitions: LigandFreq, shuffled context

**Data sources**:
- `revision_results/tool_a_ablation/tool_a_ablation_summary.csv`
- `revision_results/tool_a_ablation/tool_a_baseline_comparison.csv`
- `revision_results/multidentate/tool_a_multidentate_summary.json`
- `revision_results/multidentate/tool_a_by_denticity_final.csv`

**Tables**:
- **NEW Table S12**: Donor-field attribution ablation (full conditions from `tool_a_ablation_summary.csv`)
- **NEW Table S13**: Multidentate coverage by denticity (from `multidentate_coverage_summary.csv`)

---

### Supplementary Methods 8. Graph, GNN, and 3D equivariant baselines

**Purpose**: Full description of GNN and EGNN baselines.

**Content to include**:
1. 2D GNN architectures tested (GCN, GAT, GIN)
2. 3D EGNN architecture (equivariant graph neural network)
3. Input features, training protocol, hyperparameters
4. Task definitions: donor annotation and semantic-decoy ranking
5. Fair comparison rules (no 3D information for text baselines)
6. Leakage diagnostics: `revision_results/gnn_baselines/leakage_diagnostics.json`

**Data sources**:
- `revision_results/gnn_baselines/3d_baselines_summary.json`
- `revision_results/gnn_baselines/donor_annotation_summary.json`
- `revision_results/gnn_baselines/hard_negative_summary.json`

**Tables**:
- **NEW Table S14**: GNN architecture hyperparameters
- **NEW Table S15**: 3D EGNN donor annotation by CN (from `3d_donor_annotation_by_cn.csv`)

---

### Supplementary Methods 9. Syntax validation, repair, and semantic consistency checks

**Purpose**: Define grammar validation, repair protocol, and semantic scoring.

**Content to include**:
1. Syntax validity checker specification
2. Corruption types (existing Table S5 → renumber as S16)
3. Repair protocol: how MLM-based repair works
4. Semantic consistency scoring: CoordRep-Ranker
   - Training: contrastive learning on valid vs corrupted
   - Evaluation: AUROC on held-out decoys
5. Field ablation for semantic decoy ranking
   - Data: `revision_results/gnn_baselines/taskB_coordrep_field_ablation.csv`

**Tables**:
- **Table S16**: Corruption types (existing S5, renumbered)
- **NEW Table S17**: CoordRep-Ranker field ablation results

---

### Supplementary Methods 10. Full-CSD CoordRep-PathFinder audit

**Purpose**: Describe the 1.4M-entry CSD scan and PathFinder analysis.

**Content to include**:
1. Scan procedure and runtime
2. Waterfall filtering details (overlap with SM1 but focus on PathFinder)
3. L3 family geometry trajectory extraction
4. Boundary-record classification
5. CSD identity-layer retention statistics
   - Data: `revision_results/csd_external_summary_only/csd_identity_layer_retention.csv`
6. Transfer evaluation: ranker and Tool B on CSD
   - Data: `revision_results/csd_external_summary_only/csd_ranker_transfer_summary.csv`
   - Data: `revision_results/csd_external_summary_only/csd_toolb_transfer_summary.csv`

**Tables**:
- **NEW Table S18**: CSD identity-layer retention statistics
- **NEW Table S19**: CSD transfer evaluation summary

---

### Supplementary Methods 11. Reproducibility package and file manifest

**Purpose**: Complete reproducibility documentation.

**Content to include**:
1. Package scope and directory layout (existing SR4.1)
2. One-command reproduction instructions (existing SR4.2)
3. Deterministic seeds (existing SR4.3, Table S10 → renumber)
4. Dependencies and version pinning (existing SR4.4)
5. Data availability (existing Table S11 → renumber as S20)
6. Computational requirements (existing SR4.6)
7. **NEW**: Complete file manifest of all revision outputs

**Tables**:
- **Table S20**: Data Availability (existing S11, renumbered)
- **NEW Table S21**: File manifest for revision outputs

---

## Supplementary Results

---

### Supplementary Results 1. Dataset statistics and scope audit

**Content**:
- Full scope table: `revision_results/scope/coordrep_v1_scope_table.csv`
- Dataset composition by metal, CN, ligand count
- Split statistics verification

**Figures**:
- **Fig. S5**: Dataset composition (existing S2, renumbered)

---

### Supplementary Results 2. Canonicalization and identity robustness

**Content**:
- Empirical determinism verification results
- CoordRep-ID collision analysis across L0–L3
- CSD family identity summary: `revision_results/csd_external_summary_only/csd_family_identity_summary.csv`

**Tables**:
- **NEW Table S22**: L0–L3 collision rates

---

### Supplementary Results 3. Representation-operation challenge details

**Content**:
- Full results of the Coordination Identity Challenge (Fig. 2 expansion)
- Per-operation success/failure for each representation
- Challenge pair examples with explanations

**Data sources**:
- `revision_results/representation_gap_benchmark/coordination_identity_challenge_summary.csv`
- `revision_results/representation_gap_benchmark/challenge_pair_examples.jsonl`

**Figures**:
- **Fig. S6**: Full representation capability heatmap (expanded from Fig. 2C)

---

### Supplementary Results 4. Continuous geometry supplementary analyses

**Content**:
- CN4 SP↔Td distribution and enrichment
- CN5 Berry pathway statistics
- CN6 distortion atlas
- PathFinder boundary statistics
- L3 family geometry trajectories: `revision_results/csd_pathfinder_full/l3_family_geometry_trajectories.csv`

**Figures**:
- **Fig. S7**: Geometry region enrichment by CN (from `geometry_region_enrichment_cn*.csv`)
- **Fig. S8**: L3 family geometry trajectory examples (from `l3_family_geometry_trajectories.csv`)

---

### Supplementary Results 5. Tokenizer ablations and masked-field learning controls

**Content**:
- Factorized vs composite tokenizer full results
  - `revision_results/factorized_token/factorized_ablation_summary.csv`
  - `revision_results/factorized_token/factorized_by_task.csv`
- Donor recovery by CN: `revision_results/tool_a_ablation/tool_a_ablation_by_cn.csv`
- Donor recovery by donor element: `revision_results/tool_a_ablation/tool_a_ablation_by_donor_element.csv`
- Donor recovery by ligand frequency: `revision_results/tool_a_ablation/tool_a_ablation_by_ligand_freq.csv`

**Figures**:
- **Fig. S9**: Factorized tokenizer ablation bar chart (from Fig. 4D data)
- **Fig. S10**: Donor recovery by CN breakdown

**Tables**:
- **Table S23**: Full factorized tokenizer results by task and denticity
- **Table S24**: Donor recovery by donor element

---

### Supplementary Results 6. Donor-set recovery for ligand denticity subsets

**Content**:
- Per-denticity donor-set recovery: `revision_results/multidentate/tool_a_by_denticity_final.csv`
- Hard-error examples: `revision_results/multidentate/tool_a_multidentate_invalid_cases.csv`
- Coverage by source: `revision_results/multidentate/multidentate_by_source.csv`
- Validator failures: `revision_results/multidentate/multidentate_validator_failures.jsonl`

**Figures**:
- **Fig. S11**: Donor-set recovery by denticity (from Fig. 4E data, expanded)

**Tables**:
- **Table S25**: Multidentate donor-set recovery by denticity (full breakdown)
- **Table S26**: Hard-error examples (from `tool_a_multidentate_invalid_cases.csv`)

---

### Supplementary Results 7. Graph baseline and semantic-decoy ranking

**Content**:
- 2D GNN donor annotation full results
  - By CN: `revision_results/gnn_baselines/3d_donor_annotation_by_cn.csv`
  - By denticity: `revision_results/gnn_baselines/3d_donor_annotation_by_denticity.csv`
  - By ligand frequency: `revision_results/gnn_baselines/donor_annotation_by_ligand_freq.csv`
- 3D EGNN results: `revision_results/gnn_baselines/3d_donor_annotation_summary.csv`
- Semantic decoy ranking
  - By decoy type: `revision_results/gnn_baselines/hard_negative_by_decoy_type.csv`
  - 3D comparison: `revision_results/gnn_baselines/3d_hard_negative_by_type.csv`
  - EGNN metadata control: `revision_results/gnn_baselines/taskB_egnn_metadata_control.csv`
- CoordRep-Ranker full results: `revision_results/coordrep_ranker/ranker_summary.csv`
  - By CN: `revision_results/coordrep_ranker/ranker_by_cn.csv`
  - By boundary: `revision_results/coordrep_ranker/ranker_by_boundary.csv`
  - By decoy type: `revision_results/coordrep_ranker/ranker_by_decoy_type.csv`

**Figures**:
- **Fig. S12**: GNN donor annotation by CN/denticity
- **Fig. S13**: Semantic-decoy AUROC by decoy type (expanded Fig. 4F)

**Tables**:
- **Table S27**: Full GNN comparison table (from `gnn_comparison_summary.csv`)
- **Table S28**: EGNN metadata control ablation
- **Table S29**: CoordRep-Ranker by decoy type

---

### Supplementary Results 8. Syntax repair and semantic consistency casebook

**Content**:
- Tool B sequence baseline results: `revision_results/toolb_sequence_baselines/toolb_baseline_summary.csv`
- By corruption type: `revision_results/toolb_sequence_baselines/toolb_by_corruption_type.csv`
- CSD transfer: `revision_results/toolb_sequence_baselines/toolb_csd_transfer_summary.csv`
- Failure examples: `revision_results/toolb_sequence_baselines/toolb_failure_examples.csv`
- CSD casebook: `revision_results/csd_casebook/` (if applicable)

**Figures**:
- **Fig. S14**: Repair accuracy by corruption type (bar chart)
- **Fig. S15**: Representative correction cases (existing S5, updated)

**Tables**:
- **Table S30**: Tool B repair accuracy by corruption type
- **Table S31**: CSD transfer evaluation summary

---

### Supplementary Results 9. Full-CSD derived outputs and boundary-record statistics

**Content**:
- Full-CSD scan summary: `revision_results/csd_pathfinder_full/full_csd_pathfinder_summary.json`
- Valid records index statistics: `revision_results/csd_pathfinder_full/full_csd_valid_records_index.csv`
- PathFinder application summary: `revision_results/csd_pathfinder_full/pathfinder_application_summary.json`
- CSD identity baseline results: `revision_results/gnn_baselines/csd_identity_summary.json`
- CSD ranker transfer: `revision_results/csd_external_summary_only/csd_ranker_transfer_summary.csv`

**Figures**:
- **Fig. S16**: CSD waterfall filter (bar chart from `full_csd_filtering_waterfall.csv`)

**Tables**:
- **Table S32**: Full-CSD scan statistics (retain/reject counts)
- **Table S33**: PathFinder-derived geometry trajectory summary

---

## Supplementary Tables (consolidated numbering)

| # | Title | Source |
|---|---|---|
| S1 | Stratified split by metal series | existing |
| S2 | Stratified split by CN | existing |
| S3 | Full-CSD waterfall filter stages | `full_csd_filtering_waterfall.csv` |
| S4 | Token types | existing (was S4) |
| S5 | CoordRep-ID level definitions | new |
| S6 | Representation capability matrix | `representation_capability_matrix.csv` |
| S7 | Reference polyhedra library | existing (was S3) |
| S8 | CN5 PathFinder ridge statistics | `cn5_pathway_ridge_summary.csv` |
| S9 | Training hyperparameters | existing (was S8) |
| S10 | Training dynamics | existing (was S9) |
| S11 | Factorized tokenizer ablation | `factorized_ablation_summary.csv` |
| S12 | Donor-field attribution ablation | `tool_a_ablation_summary.csv` |
| S13 | Multidentate coverage by denticity | `multidentate_coverage_summary.csv` |
| S14 | GNN architecture hyperparameters | new |
| S15 | 3D EGNN donor annotation by CN | `3d_donor_annotation_by_cn.csv` |
| S16 | Corruption types | existing (was S5) |
| S17 | CoordRep-Ranker field ablation | `taskB_coordrep_field_ablation.csv` |
| S18 | CSD identity-layer retention | `csd_identity_layer_retention.csv` |
| S19 | CSD transfer evaluation | `csd_ranker_transfer_summary.csv` |
| S20 | Data availability | existing (was S11) |
| S21 | File manifest | new |
| S22 | L0–L3 collision rates | new (from CSD audit) |
| S23 | Factorized tokenizer by task/denticity | `factorized_by_task.csv` |
| S24 | Donor recovery by donor element | `tool_a_ablation_by_donor_element.csv` |
| S25 | Multidentate recovery by denticity | `tool_a_by_denticity_final.csv` |
| S26 | Hard-error examples | `tool_a_multidentate_invalid_cases.csv` |
| S27 | Full GNN comparison | `gnn_comparison_summary.csv` |
| S28 | EGNN metadata control ablation | `taskB_egnn_metadata_control.csv` |
| S29 | CoordRep-Ranker by decoy type | `ranker_by_decoy_type.csv` |
| S30 | Tool B repair by corruption type | `toolb_by_corruption_type.csv` |
| S31 | CSD transfer (Tool B) | `toolb_csd_transfer_summary.csv` |
| S32 | Full-CSD scan statistics | `full_csd_scan_summary.json` |
| S33 | PathFinder trajectory summary | `l3_family_geometry_trajectories.csv` |

---

## Supplementary Figures (consolidated numbering)

| # | Title | Source |
|---|---|---|
| S1 | CN4 SP↔Td atlas | `cn4_sp_td_atlas.csv` |
| S2 | CN5 Berry pathway ridge | `cn5_pathway_ridge_points.csv` |
| S3 | CN6 Oh distortion atlas | `cn6_oh_distortion_atlas.csv` |
| S4 | Pretraining dynamics | existing (was S1) |
| S5 | Dataset composition | existing (was S2) |
| S6 | Representation capability heatmap (full) | `fig_representation_gap_heatmap.csv` |
| S7 | Geometry region enrichment | `geometry_region_enrichment_cn*.csv` |
| S8 | L3 family geometry trajectories | `l3_family_geometry_trajectories.csv` |
| S9 | Factorized tokenizer ablation | Fig. 4D preview |
| S10 | Donor recovery by CN | `tool_a_ablation_by_cn.csv` |
| S11 | Donor-set recovery by denticity | Fig. 4E preview |
| S12 | GNN donor annotation by CN/denticity | `3d_donor_annotation_by_cn.csv` |
| S13 | Semantic-decoy AUROC by type | `hard_negative_by_decoy_type.csv` |
| S14 | Repair accuracy by corruption type | `toolb_by_corruption_type.csv` |
| S15 | Representative correction cases | existing (was S5) |
| S16 | CSD waterfall filter | `full_csd_filtering_waterfall.csv` |

---

## Source Data and File Manifest

List every CSV/JSON/JSONL file in `revision_results/` with:
- filename
- description
- corresponding main-text figure/table
- row/column counts

Generate from:
```bash
find revision_results/ -name "*.csv" -o -name "*.json" -o -name "*.jsonl" | sort
```

---

## Key Narrative Rules for SI

1. **Forbidden terms**: "inverse design", "ligand recommendation", "generative design", "chemical spellchecker"
2. **Preferred terms**: "masked-field donor annotation", "donor-field attribution", "syntax repair", "semantic consistency"
3. **Claims to avoid**: "CoordRep beats GNNs" — always frame as complementary
4. **3D baselines**: Acknowledge EGNN F1≈1.0 for donor annotation with 3D access; CoordRep operates without 3D
5. **Box 1 records**: Reference `box1_si_full_records.txt` for full machine-readable strings
6. **Reproducibility**: Every table must cite the source CSV/script that generated it

---

## Items to Move FROM Main Text TO SI

Based on the revision:
1. **mer-[Co(dien)(CN)₃]** record → SI (already in `box1_si_full_records.txt`)
2. Detailed tokenizer ablation breakdown (Fig. 4D expanded) → Table S11/S23
3. Per-denticity donor recovery details → Table S25
4. GNN architecture details → Table S14
5. CSD scan waterfall methodology → SM10
6. Full representation capability matrix → Table S6

---

## Priority Writing Order

1. **SM2** (add Box 1 SI records — already generated)
2. **SM4** (Coordination Identity Challenge — needs new prose)
3. **SM6–7** (masked-field + donor attribution — restructure existing)
4. **SM8** (GNN baselines — needs new section)
5. **SM10** (Full-CSD audit — needs new section)
6. **SR3** (Representation challenge details — new)
7. **SR5–7** (Tokenizer + donor + graph results — restructure existing)
8. **SR9** (CSD outputs — needs new section)
9. Remaining sections: adapt from existing with renumbering
