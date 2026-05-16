# CoordRep Representation Gap Benchmark

## Purpose

This benchmark proves that CoordRep is **not** a wrapper around CShM or CSD
curation.  It is a coordination record language that simultaneously supports:

1. **Canonical identity** – rotation/permutation-invariant identifiers
2. **Coordination stereochemistry** – explicit cis/trans, fac/mer tokens
3. **Continuous geometry** – binned CShM with raw values preserved
4. **Boundary state** – structured boundary tag, not forced one-hot shape
5. **Multi-resolution family linking** – L0–L3 hierarchy
6. **Record validation** – parseable grammar with cross-checks

## Dataset scope

| Quantity | Value |
|----------|-------|
| Total valid CoordRep records | 124,837 |
| Records with stereo constraints | 94,124 (75.4%) |
| Records with fac/mer annotation | 8,666 (6.9%) |
| Boundary records (ΔCShM < 1.0) | 17,791 (14.3%) |
| L3 groups with ≥2 distinct L1 states | 3,739 |
| L3 groups with stereo-diverse members | 745 |

## Four challenge categories

### 1. Invariance pairs (20 examples)
Same complex (polymorph pair) with identical L3 + constraint signature.
Rotation, atom-index permutation, and ligand traversal changes must NOT
change the identifier.

### 2. Stereo-different pairs (20 examples)
Same L3 ConnID but different coordination stereochemistry (cis/trans count
or fac/mer).  The representation must **separate** these.

### 3. Family geometry-state pairs (20 examples)
Same L3 ConnID but different L0/L1 geometry state.  The representation must
**link** them as the same family while **separating** their geometry states.

### 4. Boundary geometry records (20 examples)
ΔCShM < 1.0 between top-2 shape assignments.  The representation must
provide a boundary-aware field, not a forced one-hot shape label.

## Compared representations

| Representation | Description |
|----------------|-------------|
| canonical_SMILES | RDKit/CSD canonical SMILES for the full complex |
| InChI/InChIKey | IUPAC InChI (organic-focused) |
| raw_3D_coords | Cartesian coordinates from crystal structure |
| CShM_vector | Continuous Shape Measure values |
| CoordRep_full | Full CoordRep canonical string |
| CoordRep_ID (L0–L3) | Multi-resolution identity keys |

## Files

| File | Description |
|------|-------------|
| `representation_capability_matrix.csv` | 6 reps × 6 capabilities (Y/P/N + notes) |
| `coordination_identity_challenge_summary.csv` | 4 categories with pool sizes and metrics |
| `challenge_pair_examples.jsonl` | 80 concrete challenge items |
| `fig_representation_gap_heatmap.csv` | Numeric heatmap matrix (Y=1, P=0.5, N=0) |
| `README.md` | This file |

## Key result

No existing representation covers all six capabilities.
Only CoordRep simultaneously supports canonical identity, stereochemistry,
continuous geometry, boundary awareness, multi-resolution linking, and
grammar validation.

## License note

This benchmark contains CSD-derived analysis results (refcodes, CShM values,
shape labels, aggregate counts) only.  Raw coordinates are not redistributed.
