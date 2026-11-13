# CoordRep: A Canonical, Continuous, and Compositional Representation for Machine Learning in Coordination Chemistry

This `jacs-revision` branch contains the public code and evidence for the
current manuscript, **CoordRep: A Canonical, Continuous, and Compositional
Representation for Machine Learning in Coordination Chemistry**. The frozen
manuscript base is
[`release/jacs-revision-20260826`](release/jacs-revision-20260826/). The
additive second-round methods and control package is
[`release/jacs-revision-20260929`](release/jacs-revision-20260929/); it does not
replace or silently reprocess the earlier figures.
The separately released
[`pretrained weights`](release/jacs-revision-20261006-pretrained/) add model
inference assets without changing those frozen results. Availability statements
inside the older release directories describe those packages, not this later
checkpoint release.

**Property prediction: use the [November verification package](release/jacs-revision-20261111-property-verification/)
for the new checkpoint-backed results for Figure 4C,D and Tables S6–S8A.**
It covers 9 settings × 5 folds × 3 seeds. The November models were trained
anew; they are not recovered weights for the historical predictions. Original
property results and other evidence remain available in their dated directories.

**Current Figure 4 artwork and literature-referenced chemical validation:**
[`release/jacs-revision-20261114`](release/jacs-revision-20261114/) contains the
new C/D artwork (with A/B unchanged), the sanitized 84-case source register,
Supplementary Tables S16-S18 comparison data, and versioned XYZ/CIF input and
validation code. Use its Figure 4 files with the November saved-weight
statistics, not the August C/D artwork. Raw CSD files, atomic coordinate
ledgers and PyMOL scenes are not redistributed.

## Current evidence lock

- **April 2025 CSD census:** 1,371,757 entries were examined with CCDC Python
  API 3.6.0 under a frozen 96-symbol, CCDC-derived metal policy. In total,
  783,263 entries are metal-containing; 746,511 deposited structures meet the
  three-dimensional structural-audit criterion, while 36,752 remain in the
  census without a 3D structural record.
- **Record emission and source consistency:** all 746,511 targets yielded a
  schema-valid audit record, and 746,511/746,511 passed a separate reread of
  the native CSD object against the emitted source-derived signatures. This
  checks transcription under the extraction policy, not independent chemical
  correctness of the deposited assignments.
- **Structural coverage:** at least one structural metal-site record was
  retained in 741,402/746,511 entries; every in-scope metal site is structural
  in 733,004/746,511 entries. At site level, 2,574,081/2,608,448 metal sites
  are structural and 34,367 are explicitly audit-only.
- **Target classes:** 298,932 nonpolymeric single-metal, 295,164
  nonpolymeric multimetal, and 152,415 polymeric entries form a mutually
  exclusive partition of the 746,511-entry structural-audit target.
- **Typed relations:** 213,595 entries contain a shared donor group, 113,726
  contain a confirmed haptic/π site, and 129,092 contain a resolved nonzero
  lattice-translation edge. These flags are nonexclusive; 63 entries retain
  an ambiguous collective-π candidate without forcing hapticity.
- **Periodic collection:** 15,905/15,906 public CSD MOF Collection entries were
  processed. The 10,948 state-bearing entries yielded 172,332 local states;
  74,354 (43.15%) contain a nonzero translation-labelled edge. A separate
  code-only raw-CIF rerun of the pinned full archive reproduced the primary
  36,552-MID, 60,139-SID, and 23,587-split-MID atlas and all 11 primary
  public Figure 6 summary fields.

These are source-faithful structural-transcription results, not a claim that
every site has a unique exact CoordRep identity, a supported CShM vector, or a
distance-completed first coordination sphere. Donor, bridge, haptic/π, and
periodic relations in the release-wide audit are derived from the CSD-native
molecular bond graph; no distance-derived contact is added. Raw licensed CSD
coordinates and row-level ledgers are not redistributed.

The separate [literature-referenced validation](release/jacs-revision-20261114/chemical_validation/)
processes all 44 external Fe structures (47 local Fe records) and retains
all 40 selected crystal entries. Complete-reference agreement is 21/21
molecular primary donor sets (103 donors), 17/19 crystal donor sets (all
107 assigned nonmetal contacts retained, with four additional contacts at two
Cu sites), and 16/16 disorder-conditional environments (72 physical donor
images). Six directed Au-Tl relations are checked separately. The third BF4
orientation has two saved candidates but no complete independent donor
reference. Processing, invariance and chemical agreement are different
endpoints; these results are not 84/84 chemical accuracy or an expert-blind
benchmark. Source identifiers and comparison counts in the public register
are not licensed atom-resolved native CSD ledgers.

## Reproducibility boundary

The release provides source code, public aggregate/source tables, model
configurations, split and metric summaries where redistribution is permitted,
publication artwork, manifests, and checksums. These materials support
numerical auditing and rerunning the explicitly documented public analyses.
Three pretrained CoordRep checkpoints are now published with matching
tokenizers, explicit configurations and inference examples in the
[`2026-10-06 pretrained release`](release/jacs-revision-20261006-pretrained/).
They are the production MLM, the archived factorized MLM and the fine-tuned
semantic ranker used in the documented frozen inference analyses. They are
not the separately retrained Table S4 models or the Figure 4/Table S8/S8A
historical property models, whose original checkpoint files were not saved. For the
separate November property-model weights, see the verification package above. The repository does
not claim that every figure can be regenerated from aggregate CSV files
alone. Full CSD extraction requires a locally licensed April 2025 CSD installation.

The [November 14 package](release/jacs-revision-20261114/) adds executable
raw-XYZ and conditional CIF construction routes, configurable comparison
commands, and source-independent synthetic controls. Re-extraction of the
licensed crystal examples requires the reader's own authorized source data.
The publication Figure 4 can be redrawn from the public aggregate CSV and
preserved A/B artwork without CSD access or model inference.

The 2026-09-29 addendum provides an executable, code-only periodic MID/SID
canonicalizer and public-ZIP adapter; constructed L1 collision and chelate
controls; and two parameter-matched graph-only learning controls, one paired
directly with the original Table S8 base hybrid and one with the separately
trained relation-aware hybrid. In that historical run the original-base full-cohort gap/dipole MAEs
are 0.230744 eV/1.495054 D, versus 0.242297 eV/1.568577 D when its
CoordRep side fields are masked. These results are not substituted for the
distinct relation-aware protocol's absolute metrics.

The public ZIP adapter requires a locally authorized CCDC Python API and the
pinned CSD MOF Collection archive. Its full 15,906-CIF run processed 15,905
entries and reported one `site_extract:ValueError`; the secondary
corrected-CShM shape-conflict statistic was outside the minimal runner's
scope. The September graph-only package contains runnable training code and aggregate
metrics, but not the licensed or row-level inputs or trained weights needed
to reproduce the exact frozen predictions.

## Canonicalization version boundary

CoordRep **1.1.2rc3** adds ligand-local attachment-set keys, donor-attachment
orbits, exact residual-orbit enumeration, and fail-closed handling when the
metadata required to establish exchangeability is unavailable.

This change can alter molecular strings, token sequences, and L0–L3 hashes
relative to rc2. Figure 2 and Supplementary Figure S1 are the locked rc3
canonicalization evidence. Figure 6A–C report rc3 re-encoding of the exactly
unchanged frozen molecular-family cohort (3,491 families; 9,056 records); the
selection and rerun denominators are disclosed separately in the Figure 6
evidence. The historical Figure 4 and the November property verification use
pre-rc3 property inputs: their relation hybrids use rc2 inputs and their shape hybrids use
rc1 inputs, as disclosed with the figure. No Figure 4 result is relabelled as
rc3. Release-wide CSD counts are typed-record census
quantities and do not depend on final rc3 string ordering. Figure 6D–E use the
separate frozen periodic-v3 identity construction.

The XYZ input version `1.1.2rc3+xyz3.20261113` and conditional CIF input
version `cifcfg1.20261113` in the November 14 validation package have their
own source hashes and scope. They do not relabel the earlier property inputs
or alter the release-wide census or the Figure 6 atlas.

## Start here

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
python -m pip install -e "release/jacs-revision-20260826[dev]"
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full database rerun requires a licensed CCDC installation and local April
2025 CSD access. See
[`CSD_REDISTRIBUTION_NOTICE.md`](CSD_REDISTRIBUTION_NOTICE.md).
For the second-round controls, start with the addendum
[`README`](release/jacs-revision-20260929/README.md) and its
[`SHA256SUMS.txt`](release/jacs-revision-20260929/SHA256SUMS.txt).

## Repository map

| Item | Location |
|---|---|
| Frozen manuscript base | `release/jacs-revision-20260826/` |
| Second-round code and aggregate controls | `release/jacs-revision-20260929/` |
| Released pretrained checkpoints, configurations and inference | `release/jacs-revision-20261006-pretrained/` and its linked GitHub Release assets |
| All-metal CSD aggregate evidence | `release/jacs-revision-20260826/audits/full_csd/` |
| August Figure 1–6 artwork and source tables (original Figure 4C,D) | `release/jacs-revision-20260826/figures/` |
| Verification source tables for Figure 4C,D and Tables S6–S8A, model inputs and weights | [November package](release/jacs-revision-20261111-property-verification/) |
| Current Figure 4 artwork and plotting inputs | `release/jacs-revision-20261114/figures/Figure4/` |
| Literature-referenced validation and S16-S18 source comparisons | `release/jacs-revision-20261114/chemical_validation/` |
| Final validation input routes and reproduction commands | `release/jacs-revision-20261114/validation_code/` |
| Supplementary Figures S1–S2 | `release/jacs-revision-20260826/figures/Supplementary/` |
| rc3 canonicalization evidence | `release/jacs-revision-20260826/canonicalization/` |
| CoordRep 1.1.2rc3 source | `release/jacs-revision-20260826/{coordrep,coordrep_tools,brain}/` |
| Current audit-methods/Table S9 excerpt | `release/jacs-revision-20260826/supporting_information/` |

## License and citation

Code is released under the MIT License. CSD-derived outputs remain subject to
the restrictions in `CSD_REDISTRIBUTION_NOTICE.md`.
`PUBLIC_EVIDENCE_SHA256SUMS.txt` checks the current root documentation and
November 14 release indexes; each versioned package also supplies its own
file-level checksums.

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, and Compositional Representation for Machine Learning in Coordination Chemistry},
  author = {Luo, Wen-Lin and Li, Cheng-Hui},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```
