# CoordRep — JACS revision release (2026-08-19)

**Status: current and publication-facing.** This directory supersedes
`release/jacs-revision-20260806/` and synchronizes the corrected Figure 2,
Figure 3, Figure 5, and Figure 6 evidence with the revised manuscript.

## Current evidence package

| Scope | Public evidence |
|---|---|
| CoordRep implementation | `coordrep/`, `coordrep_tools/`, `brain/`, `scripts/`, `tests/` |
| Full-CSD audit implementation | `coordrep/audit/` and the `scripts/*csd_release*` utilities |
| Full-CSD aggregate evidence | `audits/full_csd/` |
| Figure 1 | `figures/Figure1/` — 3C record anatomy and scope |
| Figure 2 | `figures/Figure2/` — exact canonicalization and corrected fac/mer CShM evidence |
| Figure 3 | `figures/Figure3/` — corrected 33,863-record CShM atlas |
| Figure 4 | `figures/Figure4/` — relation/shape controls and geometric baselines |
| Figure 5 | `figures/Figure5/` — record construction and full-CSD census |
| Figure 6 | `figures/Figure6/` — L0–L3, XEYVEC, and periodic MID/SID organization |
| Frozen protocols | `protocols/` |

## Full-CSD audit lock

The April 2025 CSD contains 1,371,757 entries. The operational target is the
602,116 entries with a three-dimensional structure and at least one in-domain
d-block center.

- Schema-valid audit-record emission: **602,116/602,116**.
- Independent native-object source fidelity: **602,116/602,116**, with zero
  mismatches.
- Entries with every in-domain site structurally recorded:
  **601,402/602,116 (99.8814%)**.
- Structurally recorded sites: **1,984,062/1,986,198 (99.8925%)**.
- Entries with an in-domain shared donor group: **147,348**.
- Entries with a confirmed graph-resolved haptic/π site: **102,241**.
- Entries with a resolved nonzero lattice-translation edge: **94,452**.

The last three counts are nonexclusive. Thirty-four entries contain an
oversized collective π candidate; candidate objects do not themselves trigger
the haptic flag, although seven of those entries also contain a confirmed
haptic site.

Terminal outcomes are 293,850 resolved, 76,608 partial, 231,206
disorder-ambiguous, and 452 audit-only records. Full emission and source
fidelity therefore must not be paraphrased as 100% fully resolved structural
coverage or complete chemical-species identity.

The mutually exclusive target partition is 253,522 nonpolymeric single-center
entries without an external-metal relation, 233,622 nonpolymeric
multiple-center/external-relation entries, and 114,972 polymeric entries. Of
the polymeric entries, 114,916 have a completely resolved translation-labelled
source graph; nine retain explicit unresolved-edge evidence, and 47 are
edge-less audit-only outcomes.

## Downstream analysis tiers

The 101,878-record strict mononuclear η¹ core and the nested 84,453-record
connectivity-supported molecular tier are downstream analysis cohorts. They
are not coverage denominators and are not added to the full-CSD classes.

## Installation and checks

```bash
cd release/jacs-revision-20260819
python -m pip install -e ".[dev]"
pytest -q
cd ../..
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full database rerun additionally requires a licensed CCDC Python API and
the April 2025 CSD. Raw coordinates, structure files, and licensed row-level
ledgers are not redistributed. Aggregate summaries retain the SHA-256
commitments of the internal records and outcomes.

`PUBLIC_ARTIFACT_MANIFEST.json` and `RELEASE_SHA256SUMS.txt` cover the public
release tree. The old source-only checksum files are retired because they did
not cover figures or aggregate evidence. Checksums use canonical LF bytes for
text files and raw bytes for binary artwork, so verification is stable across
Windows checkouts and GitHub archives.
