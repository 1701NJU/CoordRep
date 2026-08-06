# CoordRep code and evidence map

The current revision is deliberately isolated under
`release/jacs-revision-20260806/`. This avoids mixing current publication
assets with the older exploratory tree.

```text
CoordRep/
├── release/jacs-revision-20260806/
│   ├── coordrep/                 # 1.1.2rc2 core, canonicalization, CShM, identity
│   ├── coordrep_tools/           # validators, graph features, baselines
│   ├── brain/                    # factorized tokenizer/ML components
│   ├── scripts/                  # reproducible experiment utilities
│   ├── tests/                    # release unit tests
│   ├── figures/Figure1–Figure6/  # vectors, public source tables, captions
│   ├── models/E3_baselines/      # SchNet, ViSNet, PaiNN-compatible controls
│   ├── protocols/                # frozen ML and periodic-v3 protocols
│   ├── audits/                   # release audit notes
│   └── RELEASE_MANIFEST.json     # locked current numbers and retired claims
├── libcoordrep/                  # historical package/results retained for provenance
├── checkpoints/                  # small public tokenizer metadata only
└── CSD_REDISTRIBUTION_NOTICE.md
```

## Where to edit or rerun

- Core record behavior: `release/jacs-revision-20260806/coordrep/`
- Multinuclear/haptic record objects: `release/jacs-revision-20260806/coordrep/v2beta/`
- Current canonicalization report: `release/jacs-revision-20260806/CANONICALIZATION_UPGRADE_REPORT.md`
- Figure source and vectors: `release/jacs-revision-20260806/figures/`
- E(3) baseline model code: `release/jacs-revision-20260806/models/E3_baselines/`
- CSD MOF periodic-v3 specification: `release/jacs-revision-20260806/protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json`

`libcoordrep/revision_results/` is retained as an audit archive. It contains
legacy Figure 5G/Rosetta and pre-v3 periodic outputs and must not be cited as
the current manuscript evidence.
