# Quickstart

## Installation

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
pip install -e .
```

## Try It Out

### 1. Generate annotated CoordRep examples
```bash
python scripts/coordrep_examples_gen.py
```

### 2. Run tokenizer audit (finds composite/factorized mismatch)
```bash
python scripts/tokenizer_audit.py
```

### 3. Multidentate ligand coverage
```bash
python scripts/multidentate_coverage.py
```

### 4. CSD external audit (requires CSD license)
```bash
python scripts/csd_external_audit.py --max_entries 1000
```

### 5. Run tests
```bash
pytest tests/ -v
```

## Key Concepts

| Concept | File |
|---------|------|
| CoordRep encoding | `coordrep/encode.py` |
| Identity keys L0–L3 | `coordrep/identity/identity_keys.py` |
| Factorized tokenizer | `brain/tokenizer.py` |
| Grammar validation | `coordrep_tools/validate.py` |
| GNN baselines | `scripts/gnn_baselines/` |

## Pre-computed Results

All revision results are available without running experiments:
```
revision_results/figure_ready/    # Figure-ready CSV tables
revision_results/rebuttal/        # Evidence matrix for rebuttal
```
