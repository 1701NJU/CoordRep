# Data Sources

## Public Datasets

| Dataset | Access | Description |
|---------|--------|-------------|
| tmQM | [Zenodo](https://zenodo.org/record/5578911) | 86k DFT-optimized TM complexes |
| COD | [cod.iucr.org](https://www.crystallography.net/cod/) | Open crystallographic data |

## CSD (Licensed)

The Cambridge Structural Database requires a CCDC license. See [CSD_LICENSE_NOTE.md](CSD_LICENSE_NOTE.md).

**Included in this repo**: aggregate statistics only (filter waterfall, rejection reasons, identity match rates).

**Not included**: raw structures, coordinates, bulk CoordRep strings, refcode-to-structure mappings.

## Pre-computed Results

All key results from revision experiments are provided as CSV/JSON:

```
revision_results/
├── figure_ready/              # Publication-ready tables
├── rebuttal/                  # Evidence matrix
├── csd_external_summary_only/ # CSD aggregate stats (no raw data)
├── gnn_baselines/             # GNN comparison results
├── tool_a_ablation/           # Tool A field attribution
├── factorized_token/          # Tokenizer ablation
├── multidentate/              # Multidentate coverage
├── toolb_sequence_baselines/  # Tool B repair results
└── csd_casebook/              # Case study reports
```

## Model Checkpoints

Not included due to size. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for expected paths and download instructions.
