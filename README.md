# CoordRep — JACS Revision Reproducibility Package

**A Canonical, Continuous, and Compositional Representation for Machine Learning in Coordination Chemistry**

Wen-Lin Luo, Cheng-Hui Li\*

State Key Laboratory of Coordination Chemistry, School of Chemistry and Chemical Engineering, Nanjing University, Nanjing 210023, P. R. China.

📧 E-mail: luowenlin862@gmail.com

---

## What This Branch Contains

This is the **`jacs-revision`** branch — the complete reproducibility package
for the JACS revision. It contains:

- **CoordRep grammar, tokenizer, serializer, validators** — the core
  representation library
- **CoordRep-ID** — four-level identity hierarchy (L0–L3)
- **Production v1 full-CSD audit** — 124,837 valid records (8.83% coverage)
- **CoordRep-v2 beta audit** — multinuclear + haptic extensions (500,299
  records, 35.40% coverage, +300.8% gain)
- **Masked-field learning diagnostics** — syntax learnability, donor
  attribution controls
- **Graph/equivariant baselines** — E-GNN, SchNet, semantic decoy ranking
- **CoordRep-Rosetta benchmark** — cross-domain MOF↔molecular retrieval
  (86.7% top-1, +18.9 pp over donor baseline)
- **Figure and table source data** — all CSV/JSON files needed to
  regenerate every manuscript figure and table

## Scope

| Level | Coverage | Description |
|---|---|---|
| **Production v1** | 124,837 records | Mononuclear, atom-resolved η¹ coordination snapshots |
| **v2 beta** | +375,462 records | Molecular multinuclear (record graphs) + haptic/π (site objects) |
| **Total** | 500,299 records | 35.40% all-CSD, 81.28% TM candidates |
| Future scope | — | Polymeric/extended networks, severe disorder, periodic CoordRep |

## Quick Navigation

| Manuscript Item | Source Data Location |
|---|---|
| Box 1 (Examples 1–3, v1) | `libcoordrep/revision_results/box1_coordrep_records/` |
| Box 1 (Examples 4–5, v2beta) | `libcoordrep/revision_results/box1_v2beta_examples/` |
| Figure 2 (CoordRep-ID) | `libcoordrep/revision_results/identity_robustness/` |
| Figure 4B–C (field learning) | `libcoordrep/revision_results/fig4_field_learning_revision/` |
| Figure 4D (graph baselines) | `libcoordrep/revision_results/gnn_baselines/` |
| Table 2 (tokenizer ablation) | `libcoordrep/revision_results/factorized_token/` |
| Table 3 (multidentate) | `libcoordrep/revision_results/multidentate/` |
| Figure 5B (scope audit) | `libcoordrep/revision_results/csd_scope_waterfall_revision/` |
| Figure 5C–D (boundary/trajectory) | `libcoordrep/revision_results/csd_pathfinder_full/` |
| Figure 5E (v2beta extension) | `libcoordrep/revision_results/full_csd_v2beta_postfix_audit/` |
| Figure 5G (Rosetta) | `libcoordrep/revision_results/coordrep_rosetta_hard_controls/` |
| Supplementary Figures | `libcoordrep/revision_results/si_figures/` |

For the complete file-by-file mapping, see **[DATA_MANIFEST.md](DATA_MANIFEST.md)**.

## Key Results

| Metric | Value |
|---|---|
| Rotation/translation invariance | 100% (4,000 trials) |
| Syntax validity after training | >99.9% |
| Tool A Top-1 accuracy | 97.8% |
| Full-CSD v1 records | 124,837 (8.83%) |
| Full-CSD v1+v2beta records | 500,299 (35.40%) |
| TM candidate coverage | 81.28% |
| Boundary geometry records | 17,791 (14.3%) |
| CN5 boundary/ridge fraction | 42.2% |
| Nontrivial L3 families | 6,358 |
| Rosetta Top-1 (hard control) | 86.7% (Δ +18.9 pp vs donor baseline) |
| Rosetta manual audit homologous | 93.3% (CoordRep) vs 56.7% (donor) |

## Reproducibility

### 1. Install environment

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
pip install -e ".[all]"
```

Requirements: Python ≥ 3.9, PyTorch ≥ 2.0, RDKit ≥ 2023.03.

### 2. Run integrity check (no GPU required)

```bash
python libcoordrep/scripts/check_revision_package.py
```

This verifies all source-data files exist, no raw CSD files are present,
and key manuscript numbers are correct.

### 3. Run unit tests

```bash
cd libcoordrep
python -m pytest tests/ -v
```

Tests cover: CoordRep parsing, roundtrip, identity keys, v2beta records,
manifest path validation.

### 4. Reproduce figures from source CSVs

All figure panels can be regenerated from precomputed CSV/JSON files in
`revision_results/` without GPU or CSD access.

### 5. Full pipeline (requires GPU + CSD license)

CSD-dependent extraction scripts require a licensed CSD installation.
See `CSD_REDISTRIBUTION_NOTICE.md`.

## CoordRep Format

A CoordRep-v1 string consists of four blocks:

```
[Metal:Fe|ox:+2|d:d6|CN:6]          ← Metal block
<ShapeBest:Oh|Class:ideal|...>       ← Shape block (continuous CShM)
{trans:L1:N:1--L2:Cl:1}             ← Stereochemistry constraints
|L1=NCCN:N:1:N:2|L2=Cl:Cl:1|       ← Ligand dictionary
```

CoordRep-v2 beta extends this with:
- **MultiMetalRecord** — metal-centered record graphs for multinuclear complexes
- **HapticRecord** — coordination-site objects for π/haptic systems

See `CODE_MAP.md` for full repository structure.

## Documentation

| File | Description |
|---|---|
| [DATA_MANIFEST.md](DATA_MANIFEST.md) | Maps every manuscript item to source files |
| [CODE_MAP.md](CODE_MAP.md) | Repository structure and script mapping |
| [CSD_REDISTRIBUTION_NOTICE.md](CSD_REDISTRIBUTION_NOTICE.md) | CSD compliance statement |

## CSD Redistribution Notice

This repository does **not** redistribute raw CSD coordinates, CIF files,
or structure files. Only permitted derived outputs are included. See
[CSD_REDISTRIBUTION_NOTICE.md](CSD_REDISTRIBUTION_NOTICE.md) for details.

## Citation

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
