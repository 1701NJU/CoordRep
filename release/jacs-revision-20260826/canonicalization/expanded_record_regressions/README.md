# Expanded-record canonicalization regressions

These programmatic regression ensembles test the typed expanded record
abstraction used for multimetal and haptic/π coordination. Fifty-five curated
multimetal templates and 33 curated haptic/π templates were each subjected to
100 nuisance variants combining container reordering and source-label
renaming, with cyclic rotation or reversal of haptic-site members where
applicable. All 5,500/5,500 multimetal and 3,300/3,300 haptic/π variants
reproduced their reference exact record.

The templates are not verified CSD structures and are not prevalence samples.
They establish invariance only for information present in the typed record
abstraction. The executable generators and exact-canonicalization tests are in
`scripts/generate_v2beta_extension.py` and `tests/test_v2beta_exact_canonicalization.py`.
