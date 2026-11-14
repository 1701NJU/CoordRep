# Fig. 4 – Masked-field probes separate ligand-intrinsic donor recovery from record-level semantic consistency

## Revision notes

1. All "inverse design" / "ligand recommendation" wording has been removed.
2. Donor masking is treated as **donor-field attribution**.
3. Composite metal/CN is **not independently evaluable** for metal and CN; marked `n.d.` in Panel D.
4. Multidentate analysis is limited to **eta1 bidentate/multidentate** chelating ligands; haptic/pi coordination is outside CoordRep v1.
5. Graph baselines are used to **separate local donor annotation from record-level semantic consistency**, not to claim CoordRep beats GNNs.
6. **No new complexes are generated.** Donor masking tests field recovery, not ligand generation.

## Panel layout (2 x 3)

| | Left | Center | Right |
|---|---|---|---|
| **Top** | A. Probe design | B. Syntax learnability | C. Donor attribution |
| **Bottom** | D. Tokenizer ablation | E. Multidentate | F. Graph baselines |

## Key numbers for caption

| Quantity | Value |
|----------|-------|
| Syntax validity at step | 2800 |
| Full-context donor Top-1 | 85.7% |
| LigandFreq lookup Top-1 | 86.2% |
| Ligand-masked donor Top-1 | 61.4% |
| Factorized metal Top-1 | 67.8% |
| Factorized CN Top-1 | 49.0% |
| Factorized joint Top-1 | 29.4% |
| Factorized donor Top-1 | 27.5% |
| Multidentate coverage | 51.2% |
| Bidentate donor Top-1 | 86.7% |
| Bidentate hard-error rate | 24.0% |
| EGNN donor F1 | 0.998 |
| CoordRep-Ranker strict AUROC | 0.800 |
| CoordRep-Ranker stereo-hard AUROC | 0.948 |

## Files

| File | Panel |
|------|-------|
| fig4A_probe_design.json | A |
| fig4B_syntax_validity_vs_steps.csv | B |
| fig4B_mask_ratio_recovery.csv | B |
| fig4B_summary.json | B |
| fig4C_donor_field_attribution.csv | C |
| fig4C_donor_field_attribution_summary.json | C |
| fig4D_tokenizer_ablation.csv | D |
| fig4D_tokenizer_ablation_summary.json | D |
| fig4E_multidentate_donor_set_recovery.csv | E |
| fig4E_multidentate_summary.json | E |
| fig4E_hard_error_examples.jsonl | E |
| fig4F_graph_baseline_donor_annotation.csv | F |
| fig4F_graph_baseline_semantic_decoys.csv | F |
| fig4F_graph_baseline_summary.json | F |
| fig4_all_caption_numbers.json | all |
| fig4_revision_readme.md | all |
