# Property-model verification · 11 November 2026

[Versioned GitHub Release: model archives, machine-readable index and SHA-256 checksums](https://github.com/1701NJU/CoordRep/releases/tag/coordrep-property-verification-20261111)

This verification concerns property prediction, not pretrained token models,
the CSD audit or the periodic analysis. The original property checkpoints were
not retained. These are **newly trained models** under the frozen protocols.
The original results are kept in the August and September evidence packages.

There are nine settings, with five graph-grouped folds and three prescribed
seeds per setting (135 selected checkpoints). Each weight archive is named
`coordrep-property-<setting>.zip`: `gine`, `base_graph_only`,
`base_shape_hidden`, `base_full`, `relation_graph_only`,
`relation_hidden`, `relation_full`, `schnet` and `visnet`.
Consult `property-weights-index.json` for the file inventory and
`SHA256SUMS-models.txt` for asset checksums. The asset inventory at the tag
defines the distributed files. No licensed CSD structures are included.

### What to download

- [Complete inputs, fixed splits and pairs, inference/scoring scripts, and reference predictions](https://github.com/1701NJU/CoordRep/releases/download/coordrep-property-verification-20261111/coordrep-property-verification-support.zip):
  the package README gives environment and extraction instructions. Inputs are
  derived from public tmQM, not from a native CSD export. GINE uses graph
  tensors checked equal for all 48,057 molecules; unused legacy side fields
  are not distributed.
- [Seven aggregate result and contrast tables, including intervals](https://github.com/1701NJU/CoordRep/releases/download/coordrep-property-verification-20261111/property-verification-aggregate-results.zip):
  cis/trans → Table S6 / Figure 4C; shape → Table S7 / Figure 4D;
  overall → Table S8; matched graph-only controls → Table S8A.
- [Weight-file index](https://github.com/1701NJU/CoordRep/releases/download/coordrep-property-verification-20261111/property-weights-index.json) and
  [weight-archive checksums](https://github.com/1701NJU/CoordRep/releases/download/coordrep-property-verification-20261111/SHA256SUMS-models.txt).

From the extracted support directory, with a model ZIP unpacked under
`weights/`, predict one complete held-out fold:

```shell
python replay.py --checkpoint weights/base_full/fold_0/seed_11 --device cuda --output members/base_full/fold_0/seed_11.csv
```

After running all 135 members as shown in the support README, use
`python score.py --members members --output scores-from-models`.
Using `score.py` without `--members` is a separate arithmetic check from
released mean predictions; it does not run any model. Full confidence intervals
require the default scoring run, not the quick diagnostics.

The support ZIP SHA-256 is
`db04c129093b7cdcea0c4016145036fc43b4b8a43db62fe045cd7b2d4dfa0525`;
the aggregate ZIP SHA-256 is
`f61fce446896200f4428f0053a02fa06609579ea3fbaaf4657f1a0e19503bac3`.
The support archive's `export_acceptance.json` reports the additional
portable-runner check: one complete fold (9,573 molecules, two properties)
for each of the nine settings. This is separate from the original acceptance
of all 135 verification weights. Re-scoring the released mean predictions
with the full resampling procedures also recovered all seven aggregate tables.

## Results from reloaded verification weights

| Comparison (mean absolute error; lower is better) | GINE | Matched field-hidden | Matched full | SchNet | ViSNet |
|---|---:|---:|---:|---:|---:|
| 34 graph-identical cis/trans pairs, dipole difference (D) | 5.869 | 5.198 | 1.340 | 0.763 | 0.530 |
| 44 graph-identical, shape-changing pairs, gap difference (eV) | 0.536 | 0.358 | 0.255 | 0.351 | 0.400 |

Restoring donor relations improved dipole-direction concordance from 73.5%
to 94.1%. Restoring shape fields improved gap-order concordance from 76.1%
to 88.6%. These are fixed, targeted subsets; the 44 pairs are selected by
the observed gap difference. Neither test establishes general model superiority.

In the full 48,057-record cohort, base-full and shape-hidden gap MAEs were
0.232512 and 0.232614 eV. The paired improvement was 0.000101 eV (95% CI
−0.000938 to 0.001146 eV; calculated before rounding). On absolute dipoles,
the base-full MAE was 1.495058 D vs 0.681919 D for SchNet and 0.587930 D
for ViSNet. The base-full and relation-full hybrids belong to distinct
protocols and should not be substituted for each other.

| Full-cohort control | All side fields masked | Full | Improvement (95% CI) |
|---|---:|---:|---:|
| Base, gap (eV) | 0.244260 | 0.232512 | 0.011747 (0.010419–0.013103) |
| Base, dipole (D) | 1.578700 | 1.495058 | 0.083641 (0.074395–0.093065) |
| Relation-aware, gap (eV) | 0.244150 | 0.226041 | 0.018109 (0.016669–0.019554) |
| Relation-aware, dipole (D) | 1.571560 | 1.337617 | 0.233944 (0.221427–0.246418) |

These four controls restore all side fields, not shape alone. Differences
and confidence intervals were computed on unrounded results.
