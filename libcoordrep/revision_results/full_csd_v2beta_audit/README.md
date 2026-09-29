# Full-CSD CoordRep-v2 Beta Audit

## Purpose

Evaluates CoordRep-v2 beta extensions on actual CSD multinuclear and
haptic/π entries to assess expanded coverage beyond v1.

## Key Principles

- **CoordRep v1 scope remains unchanged**: mononuclear atom-resolved η1.
- **v2beta results are NOT merged into v1 valid record counts.**
- **No raw CSD coordinates are exported.**
- Coverage gain is overlap-adjusted (waterfall priority ordering).
- Full support claims require manual chemical audit.

## Key Results

- Multinuclear sample: 5,000 entries,
  83.0% conversion rate
- Haptic sample: 2,000 entries,
  100.0% conversion rate
- Projected total CoordRep-compatible: 464,477
- All-CSD coverage: 32.87%
  (v1 alone: 8.83%)
- TM-candidate coverage: 75.46%
  (v1 alone: 20.3%)
- Coverage gain over v1: +272.1%
- **Claim level: B**

## Contents

| File | Description |
|---|---|
| v2beta_candidate_pool_summary.csv | Candidate pool sizes |
| coordrep_multi_full_audit.csv | Per-entry multinuclear audit |
| coordrep_haptic_full_audit.csv | Per-entry haptic audit |
| v2beta_manual_audit_sheet.csv | Manual audit sampling sheet |
| v2beta_manual_audit_summary.json | Manual audit summary |
| coordrep_expanded_coverage_summary.csv | Denominator-corrected coverage |
| v2beta_failure_taxonomy.csv | Failure mode classification |
| fig_scope_v2beta_expanded_waterfall.csv | SI figure data |
| table_v2beta_full_audit_summary.csv | SI table data |
| coordrep_v2beta_examples_for_si.md | Example records for SI |
| coordrep_v2beta_response_numbers.json | Key numbers for response |
| claim_level_recommendation.json | Go/no-go claim assessment |
| README.md | This file |

## Regeneration

```bash
cd libcoordrep
python scripts/generate_full_csd_v2beta_audit.py
```

Note: Requires CSD Python API access and takes ~2-3 hours for full scan.
