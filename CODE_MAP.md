# CoordRep code and evidence map

The authoritative revision is isolated under
`release/jacs-revision-20260826/`.

```text
CoordRep/
├── release/jacs-revision-20260826/
│   ├── coordrep/                    # rc3 grammar and attachment-aware canonicalizer
│   ├── coordrep_tools/              # validators, features, and baselines
│   ├── brain/                       # tokenizer and sequence-interface components
│   ├── canonicalization/            # Figure 2/S2 locked experiments
│   ├── scripts/                     # benchmark and audit utilities
│   ├── tests/                       # release and regression tests
│   ├── figures/Figure1–Figure6/     # current artwork, captions, source tables
│   ├── figures/Supplementary/       # Supplementary Figures S1–S2
│   ├── audits/full_csd/             # all-metal aggregate audit evidence
│   ├── supporting_information/      # current SM1 and Tables S9A–S9B excerpt
│   ├── protocols/                   # frozen ML and periodic protocols
│   └── RELEASE_MANIFEST.json        # locked claims and scope boundaries
├── libcoordrep/                     # historical package/results provenance
└── CSD_REDISTRIBUTION_NOTICE.md
```

## Where to inspect or rerun

- Molecular record behavior: `release/jacs-revision-20260826/coordrep/`
- General and legal-orbit canonicalization tests: `release/jacs-revision-20260826/canonicalization/`
- Typed CSD audit implementation: `release/jacs-revision-20260826/coordrep/audit/`
- Release-wide all-metal aggregate evidence: `release/jacs-revision-20260826/audits/full_csd/`
- Audit launch, validation, merge, and analysis utilities: `release/jacs-revision-20260826/scripts/`
- Attachment-aware repair report: `release/jacs-revision-20260826/CANONICAL_ATTACHMENT_FIX_REPORT.md`
- Figure sources and vectors: `release/jacs-revision-20260826/figures/`
- CSD MOF periodic-v3 specification: `release/jacs-revision-20260826/protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json`

`libcoordrep/revision_results/` and superseded release directories must not be
cited as current manuscript evidence.
