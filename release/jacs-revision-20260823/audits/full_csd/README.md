# April 2025 CSD full-audit evidence

This directory contains the public, aggregate evidence for the frozen
`csd-release-audit-20260818-v3` census. The database contained 1,371,757
entries; 602,116 met the three-dimensional/in-domain-d-block target rule.

## Public artifacts

- `FULL_CSD_CLAIM_READY_SUMMARY.json` — publication-ready coverage, scope,
  issue, and canonical-accounting summary.
- `SOURCE_FIDELITY_PUBLIC_SUMMARY.json` — sanitized aggregate result of the
  independent native-object reread.
- `RUN_PROVENANCE.json` — database, runtime, source-code, and internal-ledger
  hash commitments.

The row-level `RECORDS_INTERNAL.jsonl` and `ENTRY_OUTCOMES_INTERNAL.jsonl`
files are not redistributed because they are licensed CSD derivatives. Their
SHA-256 commitments are retained in the public summaries.

The accompanying implementation is under `../../coordrep/audit/` and
`../../scripts/`. A licensed CCDC installation and the frozen April 2025 CSD
are required for a full rerun.
