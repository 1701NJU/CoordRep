# CoordRep manuscript evidence and chemical validation

This package provides the Figure 4 artwork and literature-referenced chemical
validation supporting the manuscript and Supplementary Tables S16-S18. It is
an additive release on the existing `jacs-revision` branch. Original evidence
and prior output versions remain in their dated packages.

## Figure 4

Use [`figures/Figure4/`](figures/Figure4/) for the complete figure and editable
C/D panels. Panels A and B retain the submitted artwork without redrawing.
C/D use the saved-weight statistics from the
[135-model property-verification release](https://github.com/1701NJU/CoordRep/releases/tag/coordrep-property-verification-20261111).
They replace the August C/D artwork for the current property results.

| Matched comparison | Field hidden | Full | MAE reduction | Ordering concordance |
|---|---:|---:|---:|---|
| 34 cis/trans pairs, dipole difference (D) | 5.198 | 1.340 | 74.2% | 73.5% to 94.1% |
| 44 shape-changing pairs, gap difference (eV) | 0.358 | 0.255 | 28.7% | 76.1% to 88.6% |

All plotted confidence limits are the pointwise 95% attributed-graph-cluster
bootstrap endpoints from 20,000 resamples, conditional on the three-seed
ensembles. The 44-pair set is a target-conditioned descriptive comparison.
The figure does not establish universal predictive superiority.

## Chemical validation

[`chemical_validation/`](chemical_validation/) supplies the sanitized 84-case
register, qualified-reference tables, source bibliography and aggregate
outcomes. All 44 external Fe structures, 20 periodic entries and 20
disorder-labelled entries remain in the register.

| Endpoint | Result |
|---|---|
| Original external molecular inputs processed | 44/44 structures, 47 local Fe records |
| Single-Fe published CN comparison | 40/41; case 28 is documented separately |
| Complete optimized-state molecular primary donor sets | 21/21; 103/103 donors, no missing or additional donors |
| Complete crystallographic donor sets | 17/19; all 107 assigned nonmetal donor contacts retained, four additional contacts at two Cu sites |
| Directed Au-Tl neighbour relations | 6/6, separate from the nonmetal count |
| Complete disorder-conditional local environments | 16/16; 72/72 physical donor images, none missing or additional |

The third BF4 orientation has two separately retained local candidates but no
complete independent donor reference, and is not added to the 16-environment
denominator. Execution/invariance results and chemical-reference agreement
are separate endpoints; 44/44 processing is not 44/44 independently established
chemical accuracy. These literature-referenced development and regression
comparisons are not a prospective expert-blind study or an estimate of
chemical accuracy across the entire CSD.

[`validation_code/`](validation_code/) contains the versioned authors' code
and configurable input/reproduction commands. Raw-XYZ, native CSD and
conditional CIF routes have different input and validation scopes. See each
route's instructions and the separate software/source hashes. The published
algorithms and thresholds were not tuned during publication packaging.

## Files and access

`PUBLICATION_MANIFEST.json` identifies every public file, its SHA-256 and
redistribution scope; `SHA256SUMS.txt` checks the release contents. Figure
renderings are authors' derived illustrations. No PyMOL scene containing
atomic coordinates is included.

This package excludes raw CIF/RES/reflection data, native CSD structures,
third-party paper PDFs, licensed atom-level maps and vector ledgers, private
paths, access credentials and license keys. Reference DOIs and public case
identifiers permit readers to obtain sources through the original providers.
Re-extraction from CSD requires the reader's own authorized installation.
Refer to the repository's `CSD_REDISTRIBUTION_NOTICE.md` for the overall policy.

The new XYZ/CIF validation routes do not retroactively change the release-wide
CSD census, the Figure 6 atlas, or the rc1/rc2 property inputs used by the
November verification models. Existing property weights are unchanged.
