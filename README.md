# CoordRep

**A Canonical, Continuous, and Compositional Representation for Machine Learning in Coordination Chemistry**

Wen-Lin Luo, Cheng-Hui Li\*

State Key Laboratory of Coordination Chemistry, School of Chemistry and Chemical Engineering, Nanjing University, Nanjing 210023, P. R. China.

📧 E-mail: luowenlin862@gmail.com

---

## Overview

CoordRep is a canonical string representation for transition-metal coordination complexes that enables chemical language modeling. It encodes:

- **Metal identity + coordination number** as a composite token
- **Coordination geometry** as continuous shape measures (CShM), not discrete labels
- **Ligand connectivity** as canonical SMILES with donor-atom markers
- **Stereochemistry** (cis/trans, fac/mer) as explicit relational constraints

CoordRep strings are **rotation-invariant**, **permutation-invariant**, and **deterministic**: the same complex always maps to the same string regardless of coordinate frame or atom ordering.

A pretrained Masked Language Model (MLM) over CoordRep enables two downstream tools:

| Tool | Task | Description |
|------|------|-------------|
| **Tool A** | Donor Prediction | Given metal context + partial ligand field, predict missing donor atoms (design assistant) |
| **Tool B** | Structure Repair | Given a corrupted CoordRep sequence, recover valid syntax and chemical content (spellchecker) |

## Repository Structure

```
coordrep/
├── libcoordrep/
│   ├── coordrep/              # Core representation library
│   │   ├── canonical/         #   Canonicalization engine
│   │   ├── geometry/          #   CShM computation
│   │   ├── graph/             #   Molecular graph & ligand extraction
│   │   ├── io/                #   File readers (CIF, tmQM XYZ)
│   │   ├── serialize/         #   String serialization
│   │   ├── features/          #   ML feature extraction
│   │   └── validate/          #   Consistency checks
│   │
│   ├── brain/                 # Masked Language Model
│   │   ├── model.py           #   Transformer encoder
│   │   └── tokenizer.py       #   CoordRep tokenizer
│   │
│   ├── coordrep_tools/        # Downstream tools
│   │   ├── tool_a_donor.py    #   Tool A: Donor prediction
│   │   ├── tool_b_repair.py   #   Tool B: Structure repair
│   │   ├── baselines.py       #   Baseline methods
│   │   └── validate.py        #   Syntax validation
│   │
│   ├── scripts/               # Evaluation & plotting scripts
│   │
│   └── outputs/
│       └── fig5/
│           └── reproducible_data/  # CSV data for Fig. 5 reproduction
│
├── checkpoints/
│   └── pretrain_v3/           # Released model checkpoint
│
└── data/
    └── coordrep/              # Train/val/test splits (JSONL)
```

## Installation

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd coordrep
pip install -e ".[all]"
```

### Requirements

- Python ≥ 3.9
- PyTorch ≥ 2.0
- RDKit ≥ 2023.03
- NumPy, SciPy, pandas, scikit-learn, matplotlib

## Quick Start

### Convert a structure to CoordRep (command line)

```bash
# Single XYZ file
python -m libcoordrep.scripts.convert --xyz complex.xyz --bo complex.BO

# Quiet mode (output only the CoordRep string)
python -m libcoordrep.scripts.convert --xyz complex.xyz --quiet

# Batch-convert a tmQM dataset directory
python -m libcoordrep.scripts.convert --tmqm-dir tmQM-master/tmQM --output results.jsonl

# Batch-convert CIF files
python -m libcoordrep.scripts.convert --cif-dir data/cod_cif/ --output results.jsonl --limit 1000
```

### Encode a structure as CoordRep (Python API)

```python
from coordrep import encode

cc = encode(xyz_path="complex.xyz", bo_path="complex.BO")
s = cc.canonicalize().to_string()
print(s)
# [Metal:Fe|ox:+2|d:d6|CN:6]<ShapeBest:Oh|Class:ideal|Delta:3|V:0.23,8.10>{trans:L1:N:1--L2:Cl:1}|L1=NCCN|L2=Cl|
```

### Run Tool A evaluation (donor prediction)

```bash
python -m libcoordrep.scripts.run_tool_a_from_preds
```

### Run Tool B evaluation (structure repair)

```bash
python -m libcoordrep.scripts.run_tool_b_from_preds
```

### Reproduce Fig. 5 (CSV-only, no GPU required)

```bash
cd libcoordrep/outputs/fig5/reproducible_data
python plot_fig5_from_csv.py
```

## Reproducibility

We provide three levels of reproduction:

1. **Full pipeline**: Train from scratch → evaluate → generate figures (requires GPU)
2. **From checkpoint**: Load pretrained model → run Tool A/B evaluation → plot (requires GPU)
3. **CSV-only**: Regenerate all Fig. 5 panels from precomputed CSV files (CPU only, ~30 seconds)

All random seeds are fixed (see SI Section S10.3 for details).

## Key Results

| Metric | Value |
|--------|-------|
| Tool A: Overall Top-1 accuracy | 97.8% |
| Tool A: Mean ΔTop-1 vs CondFreq | +4.26 pp |
| Tool B: Validity uplift (overall) | +26.8 pp |
| Canonicalization: Rotation invariance | 100% (4,000 trials) |
| Syntax validity after training | >99.9% |

## CoordRep Format

A CoordRep string consists of four blocks:

```
[Metal:Fe|ox:+2|d:d6|CN:6]          ← Metal block
<ShapeBest:Oh|Class:ideal|...>       ← Shape block (continuous CShM)
{trans:L1:N:1--L2:Cl:1}             ← Stereochemistry constraints
|L1=NCCN:N:1:N:2|L2=Cl:Cl:1|       ← Ligand dictionary
```

## Citation

If you use CoordRep in your research, please cite:

```bibtex
@article{luo2026coordrep,
  title={CoordRep: A Canonical, Continuous, and Compositional Representation
         for Machine Learning in Coordination Chemistry},
  author={Luo, Wen-Lin and Li, Cheng-Hui},
  year={2026},
  journal={},
  note={State Key Laboratory of Coordination Chemistry,
        Nanjing University, Nanjing 210023, China}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

We thank the developers of the [tmQM dataset](https://github.com/bbskjelstad/tmqm) and the [Crystallography Open Database](https://www.crystallography.net/cod/) for making their data publicly available.
