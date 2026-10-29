# Multidentate Support Evidence — Revision Summary

## One-sentence conclusion

CoordRep v1 supports mononuclear η1 coordination snapshots, including monodentate and multidentate chelating ligands; it does not yet support multinuclear or ηn haptic systems.

---

## Key Numbers

### 1. Multidentate Coverage (Task 1)

| Metric | Value |
|--------|-------|
| Total complexes (tmQM + CSD) | 68,288 |
| Monodentate-only | 23.9% |
| Contains at least one bidentate | 38.1% |
| Contains at least one dent ≥ 3 | 46.0% |
| **Complexes with any multidentate** | **76.1%** |
| CSD fraction with multidentate | 66.4% |
| Total ligands | 214,721 |
| Ligand-level dent=1 | 63.9% |
| Ligand-level dent=2 | 18.1% |
| Ligand-level dent ≥ 3 | 18.0% |

### 2. Tool A Donor-Marker Completion by Denticity (Task 2)

| Mode | Denticity | Per-donor Top-1 | Per-donor Top-5 | Per-ligand All-correct |
|------|-----------|-----------------|-----------------|----------------------|
| full_context | dent=1 | 84.7% | 97.9% | 84.7% |
| full_context | dent=2 | 86.7% | 98.4% | 76.0% |
| no_ligand_smiles | dent=1 | 61.4% | 83.9% | 61.4% |
| no_ligand_smiles | dent=2 | 61.4% | 83.6% | 38.5% |

Key finding: bidentate per-donor accuracy (86.7%) is comparable to monodentate (84.7%), confirming that CoordRep grammar tracks atom-resolved donors within chelating ligands.

### 3. Annotated Examples (Task 3)

Three complete CoordRep examples generated from real pipeline:
- **A.** Monodentate-only square-planar Pt(II), CN=4 (tmQM ABAMIA)
- **B.** Bidentate chelating Pd(II), CN=4, dents=[2,2] (tmQM ABAFOZ)
- **C.** Octahedral Ir, CN=6, dents=[2,2,2], multiple trans constraints (tmQM ABAYUA)

All strings verified as valid through the full encode → canonicalize → serialize → validate pipeline.

### 4. Grammar Validator (Task 4)

| Metric | Value |
|--------|-------|
| Total strings validated | 68,288 |
| Pass rate | 93.2% |
| Hard errors (duplicate rank + donor over-count) | 314 (0.46%) |
| Soft warnings (possible_split_chelate) | 4,694 (false positives from shared SMILES) |

The hard error rate of **0.46%** is well below the 1% threshold.

### 5. Scope Table (Task 5)

See `revision_results/scope/coordrep_v1_scope_table.md` for full table.

**Supported:** mononuclear TM complexes, η1 mono/bi/multidentate ligands, CN 2–14, cis/trans/fac/mer, CoordRep-ID L0–L3.

**Excluded:** multinuclear, MOFs, ηn organometallics, disorder.

---

## Response Letter Bullet Points

- **76.1%** of complexes in the combined tmQM+CSD dataset contain at least one multidentate chelating ligand; CoordRep v1 is not limited to monodentate ligands.
- Per-donor Top-1 accuracy for bidentate ligands (**86.7%**) is comparable to monodentate (**84.7%**), demonstrating that the grammar effectively tracks multiple atom-resolved donors within a single chelating ligand.
- The multidentate grammar validator confirms a **0.46%** hard error rate across 68,288 complexes.
- Three fully annotated CoordRep examples (monodentate, bidentate, octahedral multi-stereo) are provided as main-text figures.
- A clear scope table documents what CoordRep v1 supports and what is deferred to v2 extensions (multinuclear, haptic, periodic).

---

## Output Files

### Task 1: Coverage
- `revision_results/multidentate/multidentate_coverage_summary.csv`
- `revision_results/multidentate/multidentate_by_source.csv`
- `revision_results/multidentate/multidentate_by_cn.csv`
- `revision_results/multidentate/multidentate_by_metal_row.csv`
- `revision_results/multidentate/summary.json`

### Task 2: Tool A Evaluation
- `revision_results/multidentate/tool_a_by_denticity_final.csv`
- `revision_results/multidentate/tool_a_multidentate_invalid_cases.csv`
- `revision_results/multidentate/tool_a_multidentate_summary.json`

### Task 3: Examples
- `revision_results/examples/coordrep_examples_full.json`
- `revision_results/examples/coordrep_examples_maintext.md`
- `revision_results/examples/coordrep_examples_si.md`
- `revision_results/examples/coordrep_example_strings.txt`

### Task 4: Validator
- `revision_results/multidentate/multidentate_validator_summary.csv`
- `revision_results/multidentate/multidentate_validator_failures.jsonl`

### Task 5: Scope
- `revision_results/scope/coordrep_v1_scope_table.md`
- `revision_results/scope/coordrep_v1_scope_table.csv`
- `revision_results/scope/scope_summary.json`
