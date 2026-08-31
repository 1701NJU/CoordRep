# CSD Scope Coverage Waterfall — CoordRep v1

## Summary

CoordRep v1 converts **124,837** CSD entries to valid records
from the **1,413,222** entries in CSD 2024.3.  The headline figure of
**8.83%** of all CSD entries may give a
misleading impression because the majority of CSD entries are purely
organic or main-group compounds with no coordination environment.

Using progressively narrower denominators:

| Denominator | Entries | Valid v1 records | Coverage |
|---|---:|---:|---:|
| All CSD entries | 1,413,222 | 124,837 | 8.83% |
| Transition-metal candidates | 615,498 | 124,837 | 20.3% |
| Mononuclear η1 coord. candidates | 217,163 | 124,837 | 57.5% |
| Intended v1 domain | 126,197 | 124,837 | 98.9% |

## Hierarchical Exclusion Waterfall

Starting from all 1,413,222 CSD 2024.3 entries, entries are removed
by their **first** (highest-priority) exclusion reason:

| Stage | Removed | Remaining |
|---|---:|---:|
| No 3D / atom-resolved structure | 78,081 | 1,335,141 |
| No transition metal | 719,643 | 615,498 |
| **→ TM candidates (denominator 2)** | — | **615,498** |
| Multinuclear / extended | 346,468 | 269,030 |
| Haptic / π (η > 1) | 51,867 | 217,163 |
| Disorder / partial occupancy | 64,098 | 153,065 |
| Incompatible bonding graph | 26,868 | — |
| Other | 0 | — |
| **→ Intended v1 domain** | — | **126,197** |
| In-domain validation failure | 1,360 | 124,837 |
| **→ Valid CoordRep v1 records** | — | **124,837** |

## Why the 8.83% is not the coverage fraction

1. **78,081** entries (5.5%) lack 3D coordinates entirely.
2. **719,643** entries (50.9%) are organic or main-group
   compounds with no transition metal.
3. Together, these account for **56.4%** of the CSD —
   they were never candidates for coordination representation.

Among the **615,498** TM candidates:
- **56.3%** are multinuclear (346,468 entries) — outside v1 scope
  but addressed by the CoordRep-Multi-v0 prototype extension.
- **8.4%** have haptic/π coordination (51,867 entries) — outside v1
  scope but addressed by the CoordRep-Haptic-v0 prototype extension.
- **10.4%** have crystallographic disorder
  (64,098 entries).

## Intended v1 conversion rate

Among the **126,197** entries within CoordRep v1's intended scope
(mononuclear, η1, CN 2–6, no disorder, no polymer), **98.9%** are
successfully converted.  The remaining **1,360**
(1.08%) fail final serialization
validation.

## Files

- `csd_scope_waterfall.csv` — Hierarchical waterfall with mutually exclusive categories
- `csd_scope_waterfall_summary.json` — Machine-readable summary with all metrics
- `csd_scope_exclusion_by_denominator.csv` — Each category as % of four denominators
- `fig5B_updated_scope_waterfall.csv` — Data for revised Figure 5B
- `csd_scope_waterfall_for_si.md` — This file
- `README.md` — Scope notes
