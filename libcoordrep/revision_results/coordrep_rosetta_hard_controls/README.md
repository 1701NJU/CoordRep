# CoordRep-Rosetta Hard Controls

## Purpose

Test whether CoordRep higher-order fields (bridge topology, haptic modes,
mu-pattern, shape) add retrieval value beyond metal/CN/donor matching.

Three sub-tasks tested:
- **A (same_metal_cn)**: candidates share same metal + CN. 41 pools.
- **B (same_metal)**: candidates share same metal. 143 pools.
- **C (same_metal_donor)**: candidates share same metal + donor set. 0 pools.

Best subtask for go/no-go: **B_same_metal** (143 pools, mean 57.1 candidates).

## Data

- **MOF nodes**: 500
- **Molecular references**: 4500 (mono=3000, multi=1500)

## Results: Best Subtask (B_same_metal)

| Method | Top-1 Correct | Top-5 Correct | AUROC | FP Rate |
|---|---|---|---|---|
| coordrep_full | 0.8671 | 0.9441 | 0.9489 | 0.2909 |
| donor_baseline | 0.6783 | 0.9301 | 0.9041 | 0.3902 |
| cshm_only | 0.3846 | 0.7343 | 0.5973 | 0.6490 |
| L3_only | 0.3636 | 0.6923 | 0.5246 | 0.7007 |
| random_same_donor | 0.1818 | 0.6014 | 0.5023 | 0.7790 |

## Results: Sub-task A (same_metal_cn, 41 pools)

| Method | Top-1 | Top-5 | AUROC |
|---|---|---|---|
| coordrep_full | 1.0000 | 1.0000 | 1.0000 |
| donor_baseline | 0.8780 | 0.9512 | 0.8358 |
| cshm_only | 0.5366 | 0.9512 | 0.5000 |
| L3_only | 0.7561 | 0.9512 | 0.5674 |
| random_same_donor | 0.4146 | 0.8537 | 0.4667 |

## Go/No-Go

| Criterion | Pass |
|---|---|
| full_beats_donor_top1 | YES |
| delta_top1_significant | YES |
| full_auroc_above_0.6 | YES |
| full_beats_random | YES |
| n_pools_sufficient | YES |

**Decision: GO_main_text**

Delta (CoordRep full - donor baseline):
- Top-1: +0.1888
- AUROC: +0.0448

## Manual Audit

| Type | N | Homologous Rate |
|---|---|---|
| coordrep_top1 | 30 | 93.3% |
| donor_top1 | 30 | 56.7% |
| hard_negative | 30 | 0.0% |

## Files

| File | Description |
|---|---|
| hard_control_summary.json | Full summary with go/no-go decision |
| hard_control_method_comparison.csv | Method comparison table |
| hard_control_detailed_results.csv | Per-pool per-method results |
| hard_pool_index.csv | Hard pool definitions |
| hard_control_manual_audit.csv | Manual audit rows |
| hard_control_manual_audit_summary.json | Audit summary |
| fig_hard_control_bars.csv | Figure-ready bar chart data |
| README.md | This file |

## No raw CSD coordinates exported.
