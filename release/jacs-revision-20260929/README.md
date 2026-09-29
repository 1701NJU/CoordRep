# CoordRep second-round evidence addendum

This directory adds focused controls to the frozen
[`2026-08-26 release`](../jacs-revision-20260826/). It does not overwrite that
release's figures, CSD counts, or property models. Its files address the
second-round questions about L1 stereochemistry, periodic MID/SID
construction, and the incremental value of CoordRep fields beyond a complete
attributed molecular graph.

| Directory | Contents | What the evidence establishes |
|---|---|---|
| [`identity/`](identity/) | Reproducible constructed cis/trans, fac/mer, and chelate controls using the frozen CoordRep 1.1.2rc3 serializer | Complete L0 retains the donor-pair map; the tested pairs share L1–L3. The controls do not estimate collision frequency in the CSD. |
| [`periodic/`](periodic/) | Code-only periodic-v3 canonicalizer, public CSD MOF ZIP adapters, worked ABADUG CIF replay, tests, and full-ZIP aggregate output | A raw-CIF replay of all 15,906 MOF Collection files: 15,905 processed, 172,332 sites, 36,552 MIDs, 60,139 SIDs, and 23,587 split MIDs. All 11 primary Figure 6 summary fields match the frozen public source. |
| [`ml/`](ml/) | Two equal-parameter graph-only ablations (the exact Table S8 base protocol and a separate relation-aware protocol), aggregate metrics, ligand-context equality and strict-context summaries | Each graph-only contrast isolates the combined CoordRep side-channel information within its own protocol. The relation-hidden contrast separately tests donor-relation fields; results from the two full-model runs are not interchangeable. |

`SHA256SUMS.txt` gives SHA-256 digests for every addendum file except itself.
The frozen base release retains its own manifests and checksums.

## Scope and redistribution

No licensed CIF, raw coordinate, CSD-derived per-site or per-refcode ledger,
record-level property table, or trained checkpoint weight is distributed
here. Public scripts can be inspected and, with authorized source inputs,
rerun. Aggregate metrics alone do not regenerate trained predictions. The
release-wide CSD source-signature comparison checks transcription under a
frozen extraction policy; it is not independent proof that every deposited
bond or periodic-image assignment is chemically correct. Structural coverage
is not synonymous with class-wide exact CoordRep canonicalization.

The periodic example and full raw-CIF runner require a locally authorized
CCDC Python API and the pinned public CSD MOF Collection ZIP. The ML runner
requires the frozen graph-grouped manifests, property table, and CoordRep
record JSONL, which are not reproduced in this code-only addendum. Follow
the per-directory READMEs before interpreting or rerunning any result.

The full periodic replay recorded one `site_extract:ValueError`. Its
minimal code path does not compute the secondary corrected-CShM
`shape_conflicting_metric_ids` field or rerun the separate 400-entry
representation-change challenge. The published aggregate contains no
licensed per-entry or per-site content.
