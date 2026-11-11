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
