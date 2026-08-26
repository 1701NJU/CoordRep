# CoordRep — JACS revision release (2026-08-26)

**Status: current and publication-facing.** This frozen directory supersedes
`release/jacs-revision-20260823/` and synchronizes the revised Figure 5 and
release-wide April 2025 CSD audit with the all-metal manuscript scope.

## Evidence map

| Scope | Public evidence |
|---|---|
| CoordRep 1.1.2rc3 implementation | `coordrep/`, `coordrep_tools/`, `brain/`, `tests/` |
| rc3 canonicalization experiments | `canonicalization/` and `CANONICAL_ATTACHMENT_FIX_REPORT.md` |
| All-metal CSD aggregates | `audits/full_csd/` |
| Current main figures | `figures/Figure1/` through `figures/Figure6/` |
| Supplementary Figures | `figures/Supplementary/` |
| Current audit methods and Tables S9A–S9B | `supporting_information/` |
| Frozen periodic protocols | `protocols/` |

## Canonicalization lock and compatibility boundary

Version 1.1.2rc3 adds source-index-free ligand attachment-set keys,
donor-specific attachment keys, donor-rank orbits, exact enumeration of legal
residual permutations, minimum-whole-record selection, and fail-closed
handling when attachment metadata is insufficient. The locked general
invariance challenge matched 300,000/300,000 variants. The separate
legal-orbit challenge matched 100,000/100,000 variants, whereas the
signature-only ablation matched 64,834/100,000.

Attachment-aware payloads can change canonical strings, token sequences, and
molecular L0–L3 hashes relative to rc2. Figure 2 and Supplementary Figure S2
are rc3 evidence. Figure 4 and the molecular-family portion of Figure 6 remain
frozen rc2 downstream evidence and are not represented as rc3 reruns.
Release-wide CSD typed-record counts do not depend on final rc3 string
ordering. Periodic MID/SID results use the separate periodic-v3 construction.

## All-metal CSD audit lock

The April 2025 CSD contains 1,371,757 entries. CCDC Python API 3.6.0 and a
frozen 96-symbol CCDC-derived metal policy identify 783,263 metal-containing
entries. The operational structural-audit target contains 746,511 deposited
structures meeting the 3D criterion; the other 36,752 metal-containing
entries remain in the census without a 3D structural record.

- Schema-valid audit records: **746,511/746,511**.
- Independent native-source fidelity: **746,511/746,511**, zero mismatches.
- Entries with at least one structural site: **741,402/746,511 (99.3156%)**.
- Entries with all metal sites structural: **733,004/746,511 (98.1906%)**.
- Structural metal sites: **2,574,081/2,608,448 (98.6825%)**; the remaining
  **34,367** sites are explicitly audit-only.
- Mutually exclusive target classes: **298,932** nonpolymeric single-metal,
  **295,164** nonpolymeric multimetal, and **152,415** polymeric entries.
- Nonexclusive complexity flags: **213,595** shared-donor entries,
  **113,726** confirmed haptic/π entries, and **129,092** entries with a
  resolved nonzero lattice-translation edge.
- Mutually exclusive terminal outcomes: **342,524** resolved, **112,519**
  partial, **286,359** disorder-ambiguous, and **5,109** all-sites audit-only.

These values establish source-faithful structural transcription. They do not
imply that every site has a unique exact CoordRep identity, a supported CShM
vector, or a chemically complete distance-derived first sphere. Donor,
bridge, haptic/π, and periodic relations use the native CSD molecular bond
graph; no distance-derived contact is added. Raw coordinates and licensed
row-level ledgers are not redistributed; their SHA-256 commitments are
retained in `audits/full_csd/`.

## Installation and checks

```bash
cd release/jacs-revision-20260826
python -m pip install -e ".[dev]"
cd ../..
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full database rerun requires licensed CCDC software and local April 2025
CSD access. `PUBLIC_ARTIFACT_MANIFEST.json` and `RELEASE_SHA256SUMS.txt` cover
the public release tree.

## Historical package

The complete 2026-08-23 package remains frozen at
`release/jacs-revision-20260823/`. It records the preceding 3D d-block audit
and is historical, not current Figure 5 or manuscript evidence.
