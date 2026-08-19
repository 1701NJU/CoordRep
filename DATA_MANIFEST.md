# Current-data manifest — JACS revision

All current publication artifacts are under
`release/jacs-revision-20260819/`. The release-level machine lock is
`RELEASE_MANIFEST.json`; `PUBLIC_ARTIFACT_MANIFEST.json` and
`RELEASE_SHA256SUMS.txt` provide file-level integrity checks.

| Manuscript item | Current public artifact | Primary evidence |
|---|---|---|
| CoordRep 1.1.2rc2 source | `release/jacs-revision-20260819/coordrep/` | molecular grammar, exact canonicalization, corrected CShM, audit records |
| Figure 1 | `figures/Figure1/` | 3C record anatomy and molecular/periodic scope |
| Figure 2 | `figures/Figure2/` | exact canonicalization, corrected EBAGAR/EBAGEV CShM values, source table |
| Figure 3 | `figures/Figure3/` | corrected 33,863-record CN4–CN6 CShM atlas and boundary fractions |
| Figure 4 | `figures/Figure4/` | matched relation/shape and E(3)-equivariant controls |
| Figure 5 | `figures/Figure5/` | chemistry-aware record construction and full-CSD audit |
| Figure 6 | `figures/Figure6/` | L0–L3 hierarchy, XEYVEC, and MID/SID organization |
| Full-CSD audit | `audits/full_csd/` | claim-ready aggregate, source-fidelity summary, provenance, and rerun instructions |
| E(3) baseline | `models/E3_baselines/` | runner/model code and compact OOF summaries |
| Periodic protocol | `protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json` | frozen canonical local-site specification |
| ML protocols | `protocols/` | corrected-CShM hybrid and ViSNet controls |

## Full-CSD numerical lock

The April 2025 audit examined 1,371,757 CSD entries and emitted 602,116
schema-valid target records. Independent source-fidelity validation passed
602,116/602,116 records. The strict entry-level structural count is
601,402/602,116, and the site-level count is 1,984,062/1,986,198. The
nonexclusive shared-donor, confirmed haptic/π, and resolved nonzero-translation
counts are 147,348, 102,241, and 94,452 entries, respectively.

The former 52,760 linked-multimetal and 38,203 haptic/π values were
pre-conversion candidate counts. They are historical and are not the current
Figure 5 or database-coverage evidence.

## Redistribution boundary

Public files contain aggregate values, protocols, source hashes, artwork, and
small synthetic/refcode test fixtures. Raw CSD coordinates, CIF/MOL/MOL2
files, and licensed entry/site ledgers are excluded. Their frozen SHA-256
commitments are retained in the public provenance summary.
