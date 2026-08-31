# CSD Scope Coverage Waterfall

## Purpose

This directory provides a **multi-denominator scope coverage analysis**
of CoordRep v1 against CSD 2024.3, responding to the reviewer question:

> "What fraction of real coordination complexes are excluded by current
> filtering criteria?"

## Key points

- **8.83%** is the fraction of *all* CSD entries converted to valid
  CoordRep v1 records.  This is **not** the fraction of real
  coordination complexes supported.
- The vast majority of excluded entries are organic / main-group (no TM)
  or lack 3D coordinates — they were never coordination candidates.
- Reviewer-facing metrics should use the **TM-candidate** or
  **intended-v1-domain** denominators.

## Scope statement

- CoordRep v1 is intentionally restricted to **mononuclear,
  atom-resolved η1 coordination snapshots**.
- Multinuclear and haptic / π systems are outside v1 scope but are
  addressed by supplementary prototype extensions
  (see `coordrep_extension_prototypes/`).

## Exclusion priority

When an entry has multiple exclusion reasons, it is classified by the
**first** (highest-priority) reason encountered in the sequential
filter pipeline:

1. No 3D / atom-resolved structure
2. No transition metal
3. Multinuclear / extended coordination
4. Haptic / π coordination (η > 1)
5. CN outside 2–6
6. Disorder / partial occupancy
7. Polymeric
8. Donor / SMILES failure
9. Other
10. In-domain CoordRep validation failure

This ensures every entry appears in exactly one exclusion category.

## No raw CSD coordinates

No raw CSD coordinates are exported.  All files contain aggregate
counts and derived statistics only.

## Regeneration

```bash
cd libcoordrep
python scripts/generate_csd_scope_waterfall.py
```
