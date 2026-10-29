# CoordRep-v2 Beta Extension

## Purpose

Demonstrates chemically principled extensions of CoordRep beyond v1:
- **CoordRep-Multi-v2beta**: metal-centered record graph for multinuclear
- **CoordRep-Haptic-v2beta**: coordination-site object for π/haptic

## Scope

- CoordRep v1 scope is **unchanged** (mononuclear, atom-resolved, η1).
- This is a **feasibility demonstration**, not full production support.
- No raw CSD coordinates are exported.

## Contents

| File | Description |
|---|---|
| coordrep_v2beta_summary.json | Machine-readable summary |
| coordrep_multi_v2beta_case_index.csv | 55 multinuclear case index |
| coordrep_haptic_v2beta_case_index.csv | 33 haptic case index |
| coordrep_multi_v2beta_records.jsonl | Full multinuclear records |
| coordrep_haptic_v2beta_records.jsonl | Full haptic records |
| coordrep_multi_v2beta_validation_report.csv | 14-gate validation |
| coordrep_haptic_v2beta_validation_report.csv | 14-gate validation |
| coordrep_v2beta_random_audit_summary.csv | Tier 2 audit projection |
| coordrep_v2beta_failure_taxonomy.csv | Failure mode categories |
| coordrep_v2beta_manual_audit_template.csv | Tier 3 audit template |
| coordrep_v2beta_examples_for_SI.md | SI narrative |
| coordrep_v2beta_methods_for_SI.md | SI methods section |
| coordrep_v2beta_response_numbers.json | Key numbers for response |
| reviewer_response_extension_paragraph.md | Response letter text |
| README.md | This file |

## Key Results

- **55** multinuclear curated cases: 100% parse+roundtrip
- **33** haptic curated cases: 100% parse+roundtrip
- Estimated coverage gain if matured: +228.1%

## Regeneration

```bash
cd libcoordrep
python scripts/generate_v2beta_extension.py
```
