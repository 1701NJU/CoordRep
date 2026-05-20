# CSD Data Redistribution Notice

## Statement

This repository **does not redistribute** raw CSD coordinates, CIF files,
structure files, or any other proprietary data from the Cambridge Structural
Database (CSD).

Only the following **permitted derived outputs** are included:

- CSD refcodes (public identifiers)
- Hashed CoordRep-ID keys (irreversible hash digests)
- Aggregate statistics (counts, rates, distributions)
- Validation labels (pass/fail, scope category)
- CShM-derived geometry descriptors (continuous shape measures)
- Scope classification labels
- Retrieval scores and ranking metrics
- Figure-source CSV tables
- CoordRep string records (canonical representation, not coordinates)

## Excluded File Types

The following file types are **not present** in this repository:

- `*.cif` — Crystallographic Information Files
- `*.res` — SHELX structure files
- `*.hkl` — Reflection data files
- `*.fcf` — Structure factor files
- Raw atomic coordinates from the CSD
- Full exported CSD structure files
- Commercial database dumps

## Reproduction of CSD-Dependent Steps

Users who wish to reproduce raw CSD-dependent extraction steps (e.g.,
`scripts/run_full_csd_scan.py`, `scripts/run_csd_pathfinder.py`) must:

1. Obtain their own licensed CSD installation
   (https://www.ccdc.cam.ac.uk/solutions/software/csd/)
2. Install the CSD Python API (`ccdc` package)
3. Run the extraction scripts against their local CSD copy

All downstream analysis (figures, tables, benchmarks) can be reproduced
from the precomputed CSV/JSON files in `revision_results/` without CSD access.

## Contact

For questions about data provenance, contact the corresponding author.
