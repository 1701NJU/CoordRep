# CoordRep Extension Prototypes

## Purpose

This directory contains **prototype extension tests** for the CoordRep
representation, specifically addressing reviewer concerns about:

1. **Multinuclear complexes** (Part A): Can CoordRep handle structures
   with 2–4 metal centers and bridging ligands?
2. **π/haptic ligands** (Part B): Can CoordRep represent coordination
   through π-fragments (η2–η6) rather than single donor atoms?

## Scope Limitations

- **CoordRep v1 remains restricted to mononuclear, atom-resolved η1
  coordination snapshots.** This is unchanged by these prototypes.
- These tests are **feasibility demonstrations only**; they do not
  claim full production support for MOFs, clusters, metallocenes,
  or all organometallics.
- The grammar extensions defined here (CoordRep-Multi-v0 and
  CoordRep-Haptic-v0) are prototypes subject to future refinement.
- **No raw CSD coordinates are exported** in any file in this directory.

## Contents

### Part A: Multinuclear (CoordRep-Multi-v0)

- `multinuclear_case_index.csv` — Curated case index
- `multinuclear_records.jsonl` — Full prototype records
- `multinuclear_validation_report.csv` — Validation checks
- `multinuclear_examples_for_si.md` — SI-ready formatted examples

### Part B: Haptic (CoordRep-Haptic-v0)

- `haptic_case_index.csv` — Curated case index
- `haptic_records.jsonl` — Full prototype records
- `haptic_validation_report.csv` — Validation checks
- `haptic_examples_for_si.md` — SI-ready formatted examples

### Combined

- `coordrep_extension_prototype_summary.json` — Machine-readable summary
- `coordrep_extension_prototypes_for_si.md` — Combined SI narrative
- `README.md` — This file

## Regeneration

```bash
cd libcoordrep
python scripts/generate_extension_prototypes.py
```

## Citation

If referencing these prototypes, cite as supplementary feasibility tests
from the CoordRep v1 manuscript. The prototypes demonstrate grammar
extensibility and do not constitute validated production tools.
