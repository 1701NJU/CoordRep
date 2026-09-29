# Reproducibility Guide

## Environment

```bash
conda env create -f environment.yml
conda activate coordrep
pip install -e .
```

Or with pip only:
```bash
pip install -r requirements.txt
pip install -e .
```

## Data Paths

| Dataset | Path | Access |
|---------|------|--------|
| tmQM | `data/tmQM/` | Public (Zenodo) |
| COD | `data/COD/` | Public |
| CSD | Requires local installation | CCDC license required |

## GPU Requirements

| Experiment | GPU | Time |
|------------|-----|------|
| Tool A ablation | 1× A100/V100 | ~30 min |
| Factorized tokenizer training | 1× GPU | ~2 hours |
| GNN baselines (donor) | 1× GPU | ~1 hour |
| GNN baselines (hard-neg) | 1× GPU | ~2 hours |
| CoordRep-Ranker training | 1× GPU | ~4 hours |
| CSD external audit | CPU only | ~1 hour |
| Identity benchmark | CPU only | ~10 min |
| Multidentate coverage | CPU only | ~5 min |

## Random Seeds

All experiments use fixed seeds for reproducibility:
- Default seed: 42
- CSD casebook: 456
- GNN training: 42
- Perturbation robustness: per-molecule indexed

## Reproducing Key Results

### 1. Identity hierarchy (no GPU needed)
```bash
python scripts/benchmark_identity_layers.py
python scripts/test_geometric_robustness.py
```

### 2. Tokenizer audit and factorized ablation
```bash
python scripts/tokenizer_audit.py
python scripts/factorized_token_ablation.py --device cuda:0
```

### 3. Tool A ablation (requires pretrained checkpoint)
```bash
python scripts/run_tool_a_ablation.py --checkpoint checkpoints/pretrain_v3/best.pt
```

### 4. Multidentate support
```bash
python scripts/multidentate_coverage.py
python scripts/multidentate_tool_a_eval.py
python scripts/multidentate_validator.py
```

### 5. GNN baselines
```bash
python scripts/gnn_baselines/run_donor_annotation.py --device cuda:0
python scripts/gnn_baselines/run_hard_negative_gnn.py --device cuda:0
```

### 6. CSD external validation (requires CSD license)
```bash
python scripts/csd_external_audit.py --max_entries 200000
python scripts/csd_identity_family_benchmark.py
python scripts/csd_toolb_transfer_eval.py
```

### 7. Tool B repair baselines
```bash
python scripts/run_toolb_baselines.py --device cuda:0
```

### 8. Figure-ready tables (no computation, just aggregation)
```bash
python scripts/make_figure_ready_tables.py
python scripts/make_rebuttal_evidence_matrix.py
```

## Model Checkpoints

Checkpoints are not included in the repository due to size. Expected paths:

| Checkpoint | Path | Size |
|------------|------|------|
| Pretrained MLM | `checkpoints/pretrain_v3/best.pt` | ~200 MB |
| Tokenizer | `checkpoints/pretrain_v3/tokenizer.json` | ~100 KB |
| CoordRep-Ranker | `checkpoints/coordrep_ranker/best_finetuned.pt` | ~200 MB |
| Factorized tokenizer | `revision_results/factorized_token/factorized_tokenizer/tokenizer.json` | ~100 KB |

If checkpoints are unavailable, all key results are provided as pre-computed CSV/JSON in `revision_results/`.

## Validation Without GPU

All figure-ready tables and evidence matrices can be inspected without any computation:
```bash
ls revision_results/figure_ready/
ls revision_results/rebuttal/
cat revision_results/figure_ready/README_figure_plan.md
```
