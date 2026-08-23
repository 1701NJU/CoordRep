# CoordRep canonical invariance stress test

Status: **PUBLIC_RC3_GENERAL_INVARIANCE_PASS**

**Implementation boundary.** The tested six-file attachment-aware repair is
the source now integrated as CoordRep 1.1.2rc3. The frozen outputs retain their
original pre-integration timestamps and hashes.

This package uses a fixed **tmQMg/PBE graph-disjoint CN4--6 benchmark cohort**.
It is not an April 2025 CSD prevalence sample. The selection joins the corrected
rc2 PBE ledger, graph/fold assignments, an independent single-transition-metal
XYZ audit, and the official individual tmQMg/PBE XYZ ZIP. Formula fallbacks
and records with nonempty rc2 validation issues are excluded. The strict
eligible denominator is 6,427 unique molecular graphs.
The locked cohort contains CN4/CN5/CN6 = 334/333/333 records, with 200 records
from each fold; `cohort_manifest.csv` records the exact fold-by-CN allocation.

## Execution

- Patch mode: `source`
- Gate: N=100, K=20 per arm; **PASS**
- Full N=1,000, K=100 per arm executed: **true**
- Locked full stress gate: **PASS**
- Arms: rigid motion only; atom permutation only; combined.
- `uncanonicalized_coordrep` is the source-order control produced from the same
  encoded object without calling `canonicalize()`.
- Raw Cartesian and ordered-distance controls retain source order and are
  quantized to six decimal places solely for exact serialization comparison.
- Full CoordRep mismatch details are written to `mismatch_variants.csv`; the
  file retains its header even when no mismatches occur.

## Results

| Arm | Representation | Matched encodings | Match fraction | Fully collapsed structures |
|---|---|---:|---:|---:|
| rigid | raw_cartesian | 0/100,000 | 0.000000 | 0/1,000 |
| rigid | ordered_distance_matrix | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| rigid | connectivity_only | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| rigid | uncanonicalized_coordrep | 99,838/100,000 | 0.998380 | 997/1,000 |
| rigid | patched_coordrep | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| permutation | raw_cartesian | 1,722/100,000 | 0.017220 | 0/1,000 |
| permutation | ordered_distance_matrix | 1,822/100,000 | 0.018220 | 0/1,000 |
| permutation | connectivity_only | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| permutation | uncanonicalized_coordrep | 18,570/100,000 | 0.185700 | 50/1,000 |
| permutation | patched_coordrep | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| combined | raw_cartesian | 0/100,000 | 0.000000 | 0/1,000 |
| combined | ordered_distance_matrix | 1,746/100,000 | 0.017460 | 0/1,000 |
| combined | connectivity_only | 100,000/100,000 | 1.000000 | 1,000/1,000 |
| combined | uncanonicalized_coordrep | 18,586/100,000 | 0.185860 | 50/1,000 |
| combined | patched_coordrep | 100,000/100,000 | 1.000000 | 1,000/1,000 |

## Auxiliary experimental fac/mer control

EBAGAR and EBAGEV have the same connectivity-only key:
**true**. Their full CoordRep
keys are distinct: **true**.
This auxiliary record-level test distinguished nuisance invariance within each
isomer from a donor-relation difference. The active Figure 2D uses the cleaner
idealized fac/mer-[Co(NH3)3Cl3] specificity control documented in the Figure 2
source table; no experimental CShM claim is made there.

## Files

- `run_canonical_stress.py`: executable harness.
- `cohort_manifest.csv`: locked 1,000-record selection and source hashes.
- `variant_summary.csv`: structure-level results for every arm and representation.
- `mismatch_variants.csv`: full-CoordRep mismatches/exceptions (header retained if empty).
- `fac_mer_specificity.csv`: auxiliary EBAGAR/EBAGEV record-level control; not
  the active Figure 2D chemical drawing.
- `summary.json`: aggregate results and gate decisions.
- `PUBLIC_PROVENANCE.json`: source hashes, transformation design, and public
  integration boundary without local filesystem paths.
- `SHA256SUMS.txt`: checksums of this evidence directory.
