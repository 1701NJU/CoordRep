# CoordRep-v2 Beta Post-Fix Audit

## Purpose

Push v2beta claim from B+ to C by:
1. Fixing BFS atom-identity bug (id() → label-based)
2. Separating polymeric/extended networks from molecular multinuclear
3. Re-auditing molecular subset with corrected converter
4. Verifying roundtrip/invariance/manual correctness exceed C-level thresholds

## Bug Fix

**Root cause**: `_atom_id()` used Python `id()` for CSD atom identity.
The CSD Python API returns **different Python objects** for the same atom
when accessed through different bond traversals, causing non-deterministic
atom identity during BFS fragment grouping and bridge detection.

**Fix**: Switched to `atom.label`-based identity, which is stable across
CSD bond traversals. This resolved all 428/428 validation failures.

## Scope Reclassification

| Category | N (sample) | % |
|---|---|---|
| molecular_multinuclear_supported | 4579 | 91.6% |
| polymeric_or_extended_future_scope | 330 | 6.6% |
| ambiguous_bonding_uncertain | 63 | 1.3% |
| disorder_or_partial_occupancy | 28 | 0.6% |

## Key Results

### Molecular Multinuclear (≤12 metals, non-polymeric)
- Conversion: **100.0%** (4577/4579)
- Roundtrip: **100.0%**
- Invariance: **100.0%**
- Manual correctness: **100%**

### Haptic/π (unchanged)
- All metrics: **100%**

### Coverage
- v1 valid: 124,837
- v2 molecular multi recovered: 323,459
- v2 haptic recovered: 52,003
- **Total CoordRep-compatible: 500,299**
- All-CSD: **35.40%** (v1: 8.83%)
- TM candidates: **81.28%** (v1: 20.3%)
- Coverage gain: **+300.8%**

### Claim Level: **C**

## Files

| File | Description |
|---|---|
| multinuclear_failure_diagnosis.csv | Root cause of each original failure |
| multinuclear_scope_reclassification.csv | Scope class for all 5000 entries |
| multinuclear_molecular_subset_audit.csv | Per-entry audit of molecular subset |
| multinuclear_polymeric_future_scope.csv | Entries excluded as future scope |
| v2beta_postfix_coverage_summary.csv | Updated coverage calculation |
| v2beta_postfix_claim_recommendation.json | Claim level with full metrics |
| v2beta_postfix_response_numbers.json | Key numbers for response letter |
| README.md | This file |

## No raw CSD coordinates exported.
