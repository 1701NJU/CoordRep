# CoordRep — JACS revision release (2026-08-23)

**Status: historical and superseded.** This directory preserves the frozen
2026-08-23 evidence package for provenance. It was superseded by
`release/jacs-revision-20260826/`, which expands the April 2025 CSD audit from
the earlier 3D d-block target to the frozen all-metal policy used in the
revised manuscript. Do not use the Figure 5 counts in this historical release
as current manuscript evidence.

## Evidence map

| Scope | Public evidence |
|---|---|
| CoordRep 1.1.2rc3 implementation | `coordrep/`, `coordrep_tools/`, `brain/`, `tests/` |
| rc3 canonicalization experiments | `canonicalization/` and `CANONICAL_ATTACHMENT_FIX_REPORT.md` |
| Full-CSD audit implementation and aggregates | `coordrep/audit/`, `scripts/`, `audits/full_csd/` |
| Current main figures | `figures/Figure1/` through `figures/Figure6/` |
| Supplementary Figures | `figures/Supplementary/` |
| Synchronized Supporting Information | `supporting_information/` |
| Frozen protocols | `protocols/` |

## Canonicalization lock and compatibility boundary

Version 1.1.2rc3 adds source-index-free ligand attachment-set keys,
donor-specific attachment keys, donor-rank orbits, exact enumeration of all
legal residual permutations, minimum-whole-record selection, and fail-closed
handling when attachment metadata is insufficient. The locked general
invariance challenge contains 1,000 graph-disjoint CN = 4–6 snapshots with 100
trials in each of three arms; CoordRep matched 300,000/300,000 variants. The
separate legal-orbit challenge matched 100,000/100,000 variants, whereas the
signature-only ablation matched 64,834/100,000.

The added attachment payload can change canonical strings, token sequences,
and molecular L0–L3 hashes relative to rc2. Figure 2 and Supplementary Figure
S2 are rc3 evidence. Figure 4 and the molecular-family portion of Figure 6
remain the frozen rc2 downstream evidence used in the manuscript; they are not
claimed as rerun under rc3. Full-CSD typed-record census quantities are
source-audit counts and do not depend on final rc3 string ordering. Periodic
MID/SID results use the separate frozen periodic-v3 construction.

## Full-CSD audit lock

The April 2025 CSD contains 1,371,757 entries. The operational target is the
602,116 entries with a three-dimensional structure and at least one atom in the
frozen set Sc–Zn, Y–Cd, La, and Hf–Hg. The frozen run used CCDC Python API
3.6.0.

- Schema-valid audit outcomes: **602,116/602,116**.
- Independent native-source fidelity: **602,116/602,116**, zero mismatches.
- Entries with all in-domain sites structurally recorded:
  **601,402/602,116 (99.8814%)**.
- Structurally recorded sites: **1,984,062/1,986,198 (99.8925%)**.
- Entries with an in-domain shared donor group: **147,348**.
- Entries with a confirmed graph-resolved haptic/π site: **102,241**.
- Entries with a resolved nonzero translation edge: **94,452**.

The last three flags are nonexclusive. Thirty-four entries contain an
oversized collective-π candidate; seven of those also contain a separate
confirmed haptic/π site. Terminal outcomes are 293,850 resolved, 76,608
partial, 231,206 disorder-ambiguous, and 452 audit-only records.

This is a predefined 3D d-block target audit, not an all-metal CSD audit.
Entries containing only other lanthanides, actinides, or main-group metal atoms
are outside the denominator. Mixed-metal target entries retain typed relations
to out-of-domain centers.

## Installation and checks

```bash
cd release/jacs-revision-20260823
python -m pip install -e ".[dev]"
pytest -q
cd ../..
python libcoordrep/scripts/check_jacs_revision_release.py
```

The database rerun requires licensed CCDC software and local April 2025 CSD
access. Raw CSD coordinates, structure files, and licensed row-level ledgers
are not redistributed. `PUBLIC_ARTIFACT_MANIFEST.json` and
`RELEASE_SHA256SUMS.txt` cover the public release tree.

## Historical artwork

The complete 2026-08-19 package remains frozen at
`release/jacs-revision-20260819/` and is not duplicated here. Its earlier
Figure 1/2/3/5/6 artifacts are superseded; the directories under `figures/`
in this release are the current manuscript versions.
