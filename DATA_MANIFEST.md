# Current-data manifest — JACS revision

The frozen manuscript figures, CSD census, and base analyses remain under
`release/jacs-revision-20260826/`. `RELEASE_MANIFEST.json`,
`PUBLIC_ARTIFACT_MANIFEST.json`, and `RELEASE_SHA256SUMS.txt` lock that base.
The additive second-round controls are under
`release/jacs-revision-20260929/`, with a separate `SHA256SUMS.txt`. Neither
release directory supersedes the other's evidence without an explicit note.
An additive pretrained-model release is documented under
`release/jacs-revision-20261006-pretrained/`; its binary checkpoints are
published as linked GitHub Release assets, with per-file and archive hashes.

| Manuscript item | Current public artifact | Primary evidence |
|---|---|---|
| CoordRep 1.1.2rc3 | `release/jacs-revision-20260826/coordrep/` | attachment-aware canonicalization and typed record grammar |
| Figure 1 | `release/jacs-revision-20260826/figures/Figure1/` | cisplatin-centered 3C record anatomy |
| Figure 2 | `release/jacs-revision-20260826/figures/Figure2/` | N = 1,000, K = 100 three-arm invariance, legal-orbit challenge, idealized fac/mer specificity |
| Figure 3 | `release/jacs-revision-20260826/figures/Figure3/` | 33,863-record CN4–CN6 CShM atlas and two CN5 examples |
| Figure 4 | `release/jacs-revision-20260826/figures/Figure4/` | frozen pre-rc3 property benchmark (relation inputs: rc2; shape inputs: rc1) and Cartesian controls |
| Figure 5 | `release/jacs-revision-20260826/figures/Figure5/` | all-metal April 2025 CSD census, structural-record audit, and typed complexity flags |
| Figure 6 | `release/jacs-revision-20260826/figures/Figure6/` | rc3 re-encoding of the unchanged frozen molecular-family cohort plus periodic-v3 MID/SID analysis |
| Supplementary Figures S1–S2 | `release/jacs-revision-20260826/figures/Supplementary/` | exact-orbit tractability (S1) and ΔS-threshold sensitivity (S2) |
| Canonicalization experiments | `release/jacs-revision-20260826/canonicalization/` | general invariance, legal-orbit screen, specificity, and regression summaries |
| Full-CSD audit | `release/jacs-revision-20260826/audits/full_csd/` | all-metal census, record audit, source fidelity, provenance, source tables, and internal-ledger commitments |
| Audit methods/Table S9 excerpt | `release/jacs-revision-20260826/supporting_information/` | synchronized SM1 and Tables S9A–S9B |
| Periodic protocol | `release/jacs-revision-20260826/protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json` | frozen canonical local-state specification |
| Pretrained sequence-model inference | `release/jacs-revision-20261006-pretrained/` plus linked Release assets | Exact archived production MLM, factorized MLM and semantic ranker parameters; tokenizer/configuration matching, source/export SHA-256 and synthetic inference checks. S5 and S14B mappings are specified per model. |

## Second-round evidence addendum

| Question | Public artifact | Evidential boundary |
|---|---|---|
| L1 cis/trans and fac/mer collisions; chelate donor pointers | `release/jacs-revision-20260929/identity/` | Constructed idealized controls run through the 1.1.2rc3 serializer; not CSD prevalence estimates. L0 differs and L1–L3 merge for each controlled stereoisomer pair. |
| Periodic MID/SID construction and full-ZIP replay | `release/jacs-revision-20260929/periodic/` | Code-only canonicalizer, ABADUG deposited-CIF example, and aggregate-only rerun of all 15,906 pinned public MOF CIFs. 15,905 processed; 172,332 sites, 36,552 MIDs, 60,139 SIDs, and 23,587 split MIDs; all 11 primary Figure 6 summary fields match. One `site_extract:ValueError` is retained. Corrected-CShM shape-conflict and the separate 400-entry re-expression challenge are outside this minimal rerun. |
| Direct Table S8 graph-only property control | `release/jacs-revision-20260929/ml/base_graph_only_*` | Equal-parameter masking of all 80–81 side fields in the exact original base-hybrid protocol. Gap MAE 0.242297→0.230744 eV and dipole MAE 1.568577→1.495054 D when fields are restored. Aggregate outputs only. |
| Separate relation-aware graph-only control | `release/jacs-revision-20260929/ml/graph_only_*` | Equal-parameter masking of all 198–200 side fields in the separately trained relation-aware protocol; aggregate full-cohort and 34-pair results. Its full model is not the Table S8 base full. |
| Ligand-context and strict-field checks | `release/jacs-revision-20260929/ml/pair_control_audit.json` and `strict_context_pair_*.csv` | Aggregate context-equality and restricted-pair summaries. No checkpoint weights or row-level CSD-derived inputs are distributed. |

The addendum intentionally excludes CIF files, coordinates, per-site or
per-refcode ledgers, record-level training tables, trained model weights, and
private-path QA artifacts. Its aggregate ML summaries do not independently
regenerate the reported out-of-fold predictions without the original input
cohort and retraining.

## Full-CSD numerical lock

The April 2025 CSD audit examined 1,371,757 entries under the frozen 96-symbol
CCDC-derived metal policy. The census contains 783,263 metal-containing
entries. Of these, 746,511 deposited structures meet the 3D structural-audit
criterion and 36,752 remain census-only.

All 746,511 target entries emitted schema-valid records and passed the
independent native-object reread. At least one structural metal site is
retained in 741,402 entries; every metal site is structural in 733,004
entries. Site accounting closes at 2,574,081 structural plus 34,367 audit-only
sites, for 2,608,448 total metal sites.

The mutually exclusive target classes are 298,932 nonpolymeric single-metal,
295,164 nonpolymeric multimetal, and 152,415 polymeric entries. Mutually
exclusive terminal outcomes are 342,524 resolved, 112,519 partial, 286,359
disorder-ambiguous, and 5,109 all-sites audit-only entries. Nonexclusive
complexity flags identify 213,595 shared-donor entries, 113,726 confirmed
haptic/π entries, and 129,092 entries with a resolved nonzero translation
edge; 63 entries retain an ambiguous collective-π candidate.

Metal-block signatures partition the target into 33,532 s-only, 54,888
p-only, 548,730 d-only, 48,805 f-only, and 60,556 mixed-block entries.

## Interpretation and redistribution boundaries

The release-wide audit establishes source-faithful typed structural records.
It does not assert an exact canonical CoordRep identity, a supported CShM
vector, or a distance-completed first sphere for every metal site. Audit
relations use the native CSD molecular bond graph, and no distance-derived
contact is added. Public files contain aggregate values, protocols, artwork,
source tables, and cryptographic commitments. Raw coordinates and licensed
row-level structure ledgers are excluded.

The public package includes three pretrained sequence-model checkpoints in
the 2026-10-06 Release assets, together with matching tokenizers,
configurations, synthetic inference examples and checksums. Every exported
tensor equals its archived source. Table S4's separately retrained
three-seed tokenizer models and Figure 4/Table S8/S8A property models did not
save checkpoint files and are not included. Code, permitted aggregate/source
tables, selected split and metric summaries, artwork and checksums support
numerical audit and the documented rerunnable analyses. The package does
not claim CSV-only regeneration of every figure; licensed row-level CSD
analyses require a local CSD installation.

The historical 2026-08-23 d-block package is retained in Git history at commit
`bba8300`; it is explicitly superseded and is not present in the active release
tree. Its Figure 5 denominators and relation counts are not current manuscript
evidence.
