# CoordRep

**CoordRep is a canonical, continuous, compositional, and multi-resolution representation layer for mononuclear η1 coordination complexes.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## v1 Scope

### Supported

- Mononuclear transition-metal complexes (all d-block metals, Sc–Hg)
- Atom-resolved η1 donors
- Monodentate and multidentate chelating ligands
- CN 2–14
- cis/trans, fac/mer, donor-relation stereo tokens
- CoordRep-ID L0–L3 multi-resolution identity hierarchy
- Continuous CShM geometry with boundary annotation

### Out of v1 Scope

- Multinuclear complexes
- MOFs / coordination polymers
- Polyoxometalates
- ηn haptic organometallics
- Positional disorder / partial occupancy
- Ambiguous weak interactions beyond threshold

See [docs/SCOPE.md](docs/SCOPE.md) for quantified exclusion statistics.

---

## Quickstart

```bash
# Install
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
pip install -e .

# Generate annotated examples
python scripts/coordrep_examples_gen.py

# Run tokenizer audit
python scripts/tokenizer_audit.py

# Multidentate coverage
python scripts/multidentate_coverage.py

# CSD external audit (requires licensed CSD installation)
python scripts/csd_external_audit.py --max_entries 1000

# Run tests
pytest tests/ -v
```

---

## Reproduce Revision Results

| Task | Script | Key Result |
|------|--------|------------|
| Identity hierarchy L0–L3 | `scripts/benchmark_identity_layers.py` | L3 = 88.9% CSD match |
| Perturbation robustness | `scripts/test_geometric_robustness.py` | σ=0.01: L0=33%, L3=94% |
| Tokenizer factorization | `scripts/factorized_token_ablation.py` | metal 13.8%→67.8% |
| Tool A field attribution | `scripts/run_tool_a_ablation.py` | full=85.7%, no_lig=61.4% |
| Multidentate support | `scripts/multidentate_coverage.py` | 76.1% multidentate |
| GNN donor baseline | `scripts/gnn_baselines/run_donor_annotation.py` | EGNN F1=0.998 |
| GNN semantic ranking | `scripts/gnn_baselines/run_hard_negative_gnn.py` | EGNN=0.50, CoordRep=0.80 |
| CSD external audit | `scripts/csd_external_audit.py` | 200k→17k (8.5%) |
| CSD identity families | `scripts/csd_identity_family_benchmark.py` | 765 families validated |
| Tool B repair baselines | `scripts/run_toolb_baselines.py` | Rule: 100% valid, MLM: 85% |
| Figure-ready tables | `scripts/make_figure_ready_tables.py` | 8 CSV/MD files |
| Rebuttal matrix | `scripts/make_rebuttal_evidence_matrix.py` | 27 comments mapped |

Pre-computed results are in `revision_results/` — no GPU needed to inspect them.

---

## CSD Note

CSD external scripts require a **licensed local CSD installation** (CCDC Python API). Raw CSD structures are **not redistributed** in this repository. Only aggregate statistics are included. See [docs/CSD_LICENSE_NOTE.md](docs/CSD_LICENSE_NOTE.md).

---

## Model Checkpoints

Checkpoints are not included in the repository due to size.

| Checkpoint | Expected Path | Size |
|------------|--------------|------|
| Pretrained MLM | `checkpoints/pretrain_v3/best.pt` | ~200 MB |
| Tokenizer | `checkpoints/pretrain_v3/tokenizer.json` | ~100 KB |
| CoordRep-Ranker | `checkpoints/coordrep_ranker/best_finetuned.pt` | ~200 MB |

All key results are available as pre-computed CSV/JSON tables in `revision_results/` even without checkpoints.

---

## Repository Structure

```
CoordRep/
├── README.md
├── LICENSE
├── CITATION.cff
├── pyproject.toml
├── environment.yml
├── requirements.txt
│
├── coordrep/                   # Core library
│   ├── encode.py               # CoordRep encoding
│   ├── canonical/              # Canonicalization
│   ├── identity/               # CoordRep-ID L0–L3
│   ├── geometry/               # CShM, shape, rel_config
│   ├── graph/                  # Graph construction
│   ├── validate/               # Grammar checks
│   └── ...
│
├── brain/                      # ML models
│   ├── tokenizer.py            # Factorized tokenizer (revised default)
│   ├── model.py                # MLM architecture
│   └── ...
│
├── coordrep_tools/             # Tools and baselines
│   ├── tokenizer_variants.py   # Legacy composite (ablation only)
│   ├── validate.py             # Grammar validator
│   ├── gnn_baselines/          # GCN, GIN, EGNN
│   ├── csd_adapter.py          # CSD API interface
│   └── ...
│
├── scripts/                    # Runnable experiment scripts
│   ├── gnn_baselines/
│   ├── make_figure_ready_tables.py
│   ├── make_rebuttal_evidence_matrix.py
│   └── ...
│
├── revision_results/           # Pre-computed results
│   ├── figure_ready/           # Figure-ready CSV tables
│   ├── rebuttal/               # Evidence matrix
│   ├── csd_external_summary_only/  # CSD stats (no raw data)
│   └── ...
│
├── docs/
│   ├── QUICKSTART.md
│   ├── DATA.md
│   ├── REPRODUCIBILITY.md
│   ├── CSD_LICENSE_NOTE.md
│   ├── REVISION_EXPERIMENTS.md
│   ├── TOKENIZER_FIX.md
│   └── SCOPE.md
│
└── tests/
    ├── test_tokenizer_factorized.py
    ├── test_identity_keys.py
    ├── test_coordrep_validator.py
    ├── test_multidentate_validator.py
    └── test_basic_encoding.py
```

---

## Citation

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, Compositional, and Multi-Resolution Representation Layer for Machine Learning in Coordination Chemistry},
  author = {Li, Cheng-Hui and Luo, Wen-Lin},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).
