# Additional graph-only and ligand-context controls

These are two additive, **within-protocol** information controls for the
frozen Figure 4 analysis. The first compares the exact Table S8 base full
hybrid (80–81 CoordRep-derived columns; 337,414–337,542 parameters) with a
new equal-parameter graph-only arm. The second compares a separately
trained relation-aware full hybrid (198–200 columns; 352,518–352,774
parameters) with its own equal-parameter graph-only arm. The base and
relation-aware full models are not the same model or run, and their
absolute errors must not be interchanged.

Each graph-only arm retains the complete atom/bond-attributed molecular
graph and its paired full model's GINE-plus-side-channel architecture, but
sets **every** train-fitted CoordRep side-field column to zero in training,
validation, and test partitions. The relation-aware protocol additionally
has a relation-hidden arm that retains its other fields and withholds the
donor-relation block. Both protocols use the frozen 48,057-record cohort,
five outer folds grouped by canonical attributed graph, three seeds
(11, 22, 33), targets, and training settings. The code-only runners are
`run_relation_aware_coordstatepairs_oof.py`, `run_coordstatepairs_oof.py`,
and `run_downstream_graph_benchmark.py`.

The **direct Table S8 base-model control** is:

| Quantity | Base graph-only | Original Table S8 full | Reduction (95% graph-cluster CI) |
|---|---:|---:|---:|
| Full-cohort HOMO–LUMO gap MAE | 0.242297 eV | 0.230744 eV | 0.011553 eV (0.010243–0.012876) |
| Full-cohort dipole MAE | 1.568577 D | 1.495054 D | 0.073524 D (0.064293–0.082738) |
| 34 graph-identical cis/trans pair-difference MAE | 5.868915 D | 5.258340 D | 0.610575 D (0.264760–1.036200) |

The **separate relation-aware control** is:

| Quantity | Graph-only | Full relation-aware | Difference |
|---|---:|---:|---:|
| Full-cohort HOMO–LUMO gap MAE | 0.244132 eV | 0.227769 eV | 0.016363 eV lower |
| Full-cohort dipole MAE | 1.579948 D | 1.340789 D | 0.239159 D lower |
| 34 graph-identical cis/trans pair-difference MAE | 5.868915 D | 1.327899 D | 4.541016 D lower |

The 34 pairs share the attributed graph and all recorded ligand payload,
denticity, and donor-element fields within each pair. Their graph-only
predictions tie by construction. Relation-hidden versus relation-aware
pair-difference MAE is 5.218489 versus 1.327899 D. The coarse relation
fingerprint distinguishes 31 of the 34 pairs; three are nonidentifiable
from that fingerprint alone. The frozen spin annotation is uninformative,
so these checks do not claim physical spin matching.

The 44 shape-changing pairs are target-conditioned on the observed gap
difference. Only 13 of them also match every recorded nongeometry parser
field. In that restricted post-hoc subset, the earlier shape-hidden minus
shape-aware MAE reduction is 0.238915 eV, with a graph-cluster bootstrap
95% interval of −0.062019 to 0.381499 eV. This is imprecise and is **not**
an independent positive-effect finding. The 110-pair unconditioned set
has 24 strictly matched nongeometry-context pairs; its corresponding
reduction is 0.111119 eV (interval −0.055742 to 0.281857 eV).

The `base_graph_only_*` and `graph_only_*` metrics and audit files,
`pair_control_audit.json`, and the two `strict_context_pair_*.csv` files
are aggregate-only artifacts. They do
not contain ligand strings, record IDs, CSD coordinates, row-level
predictions, or trained weights. The comparison audit records the exact
source hashes and confirms all 15 fold×seed members and equal parameter
counts **within each protocol**. The printed intervals are descriptive graph-cluster bootstraps
conditional on the frozen ensembles and are not multiple-test adjusted.

## Reuse with authorized inputs

The training scripts are provided for inspection and rerunning with a
user-supplied, authorized cohort. The exact frozen property table,
record-level CoordRep JSONL, fold/pair manifests, and trained weights are
**not** bundled here. Consequently, this aggregate package does not by
itself reproduce the reported predictions. From the repository root, with
PyTorch, PyG, RDKit, pandas, scikit-learn, and NumPy installed, the matched
relation-aware graph-only arm can be invoked as follows after supplying
those inputs:

```powershell
python release/jacs-revision-20260929/ml/run_relation_aware_coordstatepairs_oof.py `
  --record-jsonl 'C:\path\to\authorized-records.jsonl' `
  --property-csv 'C:\path\to\authorized-properties.csv' `
  --fold-manifest 'C:\path\to\frozen-record-folds.csv' `
  --pair-manifest 'C:\path\to\frozen-pairs.csv' `
  --models hybrid_mask_all --field-mode 'hybrid_mask_all=*' `
  --out 'C:\path\to\new-output-directory' --epochs 40 --patience 8 `
  --seeds 11 22 33
```

The same wrapper takes corresponding options for its full and
relation-hidden arms. The base-model graph-only arm uses
`run_coordstatepairs_oof.py` with its frozen base inputs and
`--models hybrid_mask_all --field-mode 'hybrid_mask_all=*'`. No checkpoint
files are present in this addendum.
