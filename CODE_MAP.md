# CoordRep code and evidence map

The authoritative revision is isolated under
`release/jacs-revision-20260819/`. Earlier exploratory results remain available
only as provenance and are not current manuscript evidence.

```text
CoordRep/
├── release/jacs-revision-20260819/
│   ├── coordrep/                 # core grammar, CShM, identity, audit records
│   ├── coordrep_tools/           # validators, features, and baselines
│   ├── brain/                    # tokenizer and sequence-interface components
│   ├── scripts/                  # figure, benchmark, and full-CSD audit utilities
│   ├── tests/                    # release and audit regression tests
│   ├── figures/Figure1–Figure6/  # vectors, public source tables, captions
│   ├── audits/full_csd/           # aggregate evidence and run provenance
│   ├── models/E3_baselines/      # SchNet, ViSNet, and compatible controls
│   ├── protocols/                # frozen ML and periodic protocols
│   └── RELEASE_MANIFEST.json     # locked current claims and boundaries
├── release/jacs-revision-20260806/ # superseded entry point only
├── libcoordrep/                  # historical package/results provenance
└── CSD_REDISTRIBUTION_NOTICE.md
```

## Where to edit or rerun

- Molecular record behavior: `release/jacs-revision-20260819/coordrep/`
- Full-CSD typed audit records: `release/jacs-revision-20260819/coordrep/audit/`
- Full-CSD launch, validation, merge, and analysis: `release/jacs-revision-20260819/scripts/`
- Current canonicalization report: `release/jacs-revision-20260819/CANONICALIZATION_UPGRADE_REPORT.md`
- Figure sources and vectors: `release/jacs-revision-20260819/figures/`
- CSD MOF periodic-v3 specification: `release/jacs-revision-20260819/protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json`

`libcoordrep/revision_results/` and git history retain superseded assets. They
must not be cited as evidence for the current manuscript.
