# CoordRep code and evidence map

The authoritative revision is isolated under
`release/jacs-revision-20260823/`. Earlier frozen releases and exploratory
results remain only as provenance.

```text
CoordRep/
├── release/jacs-revision-20260823/
│   ├── coordrep/                    # rc3 grammar and attachment-aware canonicalizer
│   ├── coordrep_tools/              # validators, features, and baselines
│   ├── brain/                       # tokenizer and sequence-interface components
│   ├── canonicalization/            # Figure 2/S2 locked experiments
│   ├── scripts/                     # benchmark and full-CSD audit utilities
│   ├── tests/                       # release and regression tests
│   ├── figures/Figure1–Figure6/     # current artwork, captions, source tables
│   ├── figures/Supplementary/       # Supplementary Figures S1–S2
│   ├── audits/full_csd/             # aggregate audit evidence and provenance
│   ├── supporting_information/      # synchronized SI artifacts
│   ├── protocols/                   # frozen ML and periodic protocols
│   └── RELEASE_MANIFEST.json        # locked claims and compatibility boundaries
├── release/jacs-revision-20260819/  # superseded frozen rc2 release
├── release/jacs-revision-20260806/  # superseded entry point
├── libcoordrep/                     # historical package/results provenance
└── CSD_REDISTRIBUTION_NOTICE.md
```

## Where to edit or rerun

- Molecular record behavior: `release/jacs-revision-20260823/coordrep/`
- General and legal-orbit canonicalization tests: `release/jacs-revision-20260823/canonicalization/`
- Full-CSD typed audit records: `release/jacs-revision-20260823/coordrep/audit/`
- Full-CSD launch, validation, merge, and analysis: `release/jacs-revision-20260823/scripts/`
- Attachment-aware repair report: `release/jacs-revision-20260823/CANONICAL_ATTACHMENT_FIX_REPORT.md`
- Figure sources and vectors: `release/jacs-revision-20260823/figures/`
- CSD MOF periodic-v3 specification: `release/jacs-revision-20260823/protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json`

`libcoordrep/revision_results/` and the superseded release directories must not
be cited as current manuscript evidence.
