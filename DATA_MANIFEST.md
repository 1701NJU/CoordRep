# Current-data manifest — JACS revision

All current publication assets are under
`release/jacs-revision-20260806/`. The table below is the short reviewer-facing
map; the per-figure directories contain the exact vector files, source tables,
and captions.

| Manuscript item | Current public artifact | Primary evidence |
|---|---|---|
| CoordRep 1.1.2rc2 source | `release/jacs-revision-20260806/coordrep/` | canonical and v2-beta tests |
| Figure 1 | `figures/Figure1/` | 3C record anatomy; CoSyMLib references; molecular/periodic scope |
| Figure 2 | `figures/Figure2/` | exact canonicalization and identity source files |
| Figure 3 | `figures/Figure3/` | CN4/CN5/CN6 CShM atlas tables and preview |
| Figure 4 | `figures/Figure4/` | relation/shape source CSV, SchNet/ViSNet controls, vectors |
| Figure 5 | `figures/Figure5/` | current CSD scope source CSV and final vector artwork |
| Figure 6 | `figures/Figure6/` | periodic-v3 source CSV, caption, final vector artwork |
| E(3) baseline | `models/E3_baselines/` | runner/model code and compact OOF summaries |
| Periodic protocol | `protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json` | frozen canonical local-site specification |
| ML protocol | `protocols/CORRECTED_CSHM_HYBRID_OOF_PROTOCOL_v1.json` and `protocols/VISNET_COORDSTATEPAIRS_PROTOCOL_v1.json` | frozen folds, targets, and controls |

## Redistribution boundary

The public source tables contain aggregate values and, where needed for audit,
CSD refcodes or hashes. Raw CSD coordinates, CIF/MOL/MOL2 files, commercial
database exports, and internal-licensed entry/site tables are excluded. The
full local CSD audit can be rerun only with a licensed CCDC installation.

## Current numerical lock

The authoritative machine-readable lock is
`release/jacs-revision-20260806/RELEASE_MANIFEST.json`. In particular, use the
periodic-v3 values 15,905 processed entries, 10,948 state-bearing entries,
172,332 local states, 74,354 translated-edge states (43.15%), 36,552 metric
IDs, and 60,139 stereo IDs. The earlier 15,901/172,783/10,957/43.11% package is
retired.
