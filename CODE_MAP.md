# CoordRep Code Map

## Repository Layout

```
CoordRep/
├── README.md                          ← Reviewer-facing overview
├── DATA_MANIFEST.md                   ← Maps every manuscript item to source files
├── CSD_REDISTRIBUTION_NOTICE.md       ← CSD compliance statement
├── CODE_MAP.md                        ← This file
├── LICENSE
├── pyproject.toml
│
├── checkpoints/
│   └── pretrain_v3/                   ← Released MLM checkpoint (tokenizer + config)
│
└── libcoordrep/                       ← Main package
    ├── environment.yml
    ├── requirements.txt
    ├── pyproject.toml
    │
    ├── coordrep/                      ← Core representation library
    │   ├── core.py                    ← CoordRep data structures
    │   ├── encode.py                  ← XYZ/CIF → CoordRep conversion
    │   ├── canonical/                 ← Canonicalization engine
    │   ├── geometry/                  ← CShM (continuous shape measures)
    │   ├── graph/                     ← Molecular graph & ligand extraction
    │   ├── identity/                  ← CoordRep-ID hierarchy (L0–L3)
    │   ├── io/                        ← File readers (CIF, XYZ, tmQM)
    │   ├── serialize/                 ← String serialization
    │   ├── features/                  ← ML feature extraction
    │   ├── validate/                  ← Consistency validators
    │   └── v2beta/                    ← v2-beta extensions
    │       ├── core.py                ←   MultiMetalRecord + HapticRecord
    │       ├── csd_v2beta_adapter.py  ←   CSD → v2beta conversion
    │       ├── serialize.py           ←   v2beta serialization
    │       ├── validate.py            ←   14-gate validation
    │       └── canonicalize.py        ←   v2beta canonicalization
    │
    ├── brain/                         ← Masked Language Model
    │   ├── model.py                   ← Transformer encoder
    │   └── tokenizer.py               ← CoordRep tokenizer (factorized)
    │
    ├── coordrep_tools/                ← Downstream tools
    │   ├── tool_a_donor.py            ← Tool A: donor prediction
    │   ├── tool_b_repair.py           ← Tool B: structure repair
    │   ├── baselines.py               ← Baseline methods
    │   └── validate.py                ← Syntax validation
    │
    ├── scripts/                       ← All generation & evaluation scripts
    │   ├── export_box1_coordrep_records.py        ← Box 1
    │   ├── test_geometric_robustness.py           ← Figure 2
    │   ├── benchmark_identity_layers.py           ← Figure 2
    │   ├── csd_identity_family_benchmark.py       ← Figure 2
    │   ├── export_fig4_field_learning.py           ← Figure 4
    │   ├── factorized_token_ablation.py           ← Table 2
    │   ├── multidentate_tool_a_eval.py            ← Table 3
    │   ├── run_full_csd_scan.py                   ← Figure 5B
    │   ├── generate_csd_scope_waterfall.py        ← Figure 5B
    │   ├── run_csd_pathfinder.py                  ← Figure 5C–D
    │   ├── export_fig5c_case_xyz.py               ← Figure 5C
    │   ├── export_fig5d_family_trajectory.py      ← Figure 5D
    │   ├── build_merged_fig5.py                   ← Figure 5 merged
    │   ├── generate_full_csd_v2beta_audit.py      ← Figure 5E (v2beta)
    │   ├── generate_v2beta_postfix_audit.py       ← Figure 5E (postfix)
    │   ├── generate_rosetta_homology.py           ← Figure 5G (Rosetta)
    │   ├── generate_rosetta_hard_controls.py      ← Figure 5G (hard ctrl)
    │   ├── gnn_baselines/                         ← Figure 4D baselines
    │   ├── si_figures/                            ← Supplementary figures
    │   └── check_revision_package.py              ← Package integrity test
    │
    ├── tests/                         ← Unit tests
    │   ├── test_basic_encoding.py
    │   ├── test_coordrep_validator.py
    │   ├── test_identity_keys.py
    │   ├── test_invariance.py
    │   ├── test_multidentate_validator.py
    │   ├── test_tokenizer_factorized.py
    │   ├── test_v2beta_records.py     ← v2beta validation tests
    │   └── test_manifest_paths.py     ← Manifest path checker
    │
    ├── revision_results/              ← All revision source data
    │   ├── box1_coordrep_records/     ← Box 1 (Examples 1–3, v1)
    │   ├── box1_v2beta_examples/      ← Box 1 (Examples 4–5, v2beta)
    │   ├── identity_robustness/       ← Figure 2
    │   ├── factorized_token/          ← Table 2
    │   ├── fig4_field_learning_revision/  ← Figure 4
    │   ├── gnn_baselines/             ← Figure 4D
    │   ├── tool_a_ablation/           ← Figure 4C controls
    │   ├── multidentate/              ← Table 3
    │   ├── csd_scope_waterfall_revision/  ← Figure 5B
    │   ├── csd_pathfinder_full/       ← Figure 5C–D
    │   ├── full_csd_v2beta_audit/     ← Figure 5E (initial)
    │   ├── full_csd_v2beta_postfix_audit/ ← Figure 5E (postfix)
    │   ├── coordrep_v2_beta_extension/    ← v2beta curated cases
    │   ├── coordrep_rosetta_active_site_homology/ ← Figure 5G pilot
    │   ├── coordrep_rosetta_hard_controls/        ← Figure 5G hard ctrl
    │   ├── figure_ready_merged_fig5/  ← Figure 5 merged panels
    │   ├── si_figures/                ← Supplementary figure sources
    │   ├── figure_source_data/        ← All figure CSV source data
    │   └── manuscript_tables/         ← All table CSV source data
    │
    └── docs/                          ← Documentation
        ├── grammar_specification.md
        ├── coordrep_id.md
        └── v2beta_scope.md
```

## Manuscript → Code/Data Mapping

| Manuscript Item | Primary Script | Output Directory |
|---|---|---|
| Box 1 (Examples 1–3) | `export_box1_coordrep_records.py` | `revision_results/box1_coordrep_records/` |
| Box 1 (Examples 4–5) | `generate_v2beta_postfix_audit.py` | `revision_results/box1_v2beta_examples/` |
| Figure 2 | `test_geometric_robustness.py`, `csd_identity_family_benchmark.py` | `revision_results/identity_robustness/` |
| Figure 4B | `export_fig4_field_learning.py` | `revision_results/fig4_field_learning_revision/` |
| Figure 4C | `export_fig4_field_learning.py` | `revision_results/fig4_field_learning_revision/` |
| Figure 4D | `gnn_baselines/` scripts | `revision_results/gnn_baselines/` |
| Table 2 | `factorized_token_ablation.py` | `revision_results/factorized_token/` |
| Table 3 | `multidentate_tool_a_eval.py` | `revision_results/multidentate/` |
| Figure 5B | `generate_csd_scope_waterfall.py` | `revision_results/csd_scope_waterfall_revision/` |
| Figure 5C | `run_csd_pathfinder.py` | `revision_results/csd_pathfinder_full/` |
| Figure 5D | `export_fig5d_family_trajectory.py` | `revision_results/csd_pathfinder_full/` |
| Figure 5E | `generate_v2beta_postfix_audit.py` | `revision_results/full_csd_v2beta_postfix_audit/` |
| Figure 5G | `generate_rosetta_hard_controls.py` | `revision_results/coordrep_rosetta_hard_controls/` |
| Supp. Figures | `scripts/si_figures/` | `revision_results/si_figures/` |
