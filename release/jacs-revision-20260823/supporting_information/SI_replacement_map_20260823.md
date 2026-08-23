# SI replacement map — 2026-08-23

Source preserved unchanged: `Supplementary Information - clean.docx`
New synchronized file: `Supplementary Information - synchronized 20260823.docx`

## Global metadata and audit scope

- Replaced the title with **CoordRep: A Canonical, Continuous, and
  Compositional Representation for Coordination Chemistry**.
- Corrected the frozen release-wide CSD runtime from CCDC Python API 3.4.0 to
  **3.6.0**.
- Defined the frozen in-domain element set explicitly as **Sc–Zn, Y–Cd, La,
  and Hf–Hg**.
- Stated that entries containing only Ce–Lu, actinides, or main-group metal
  atoms are outside the 602,116-entry denominator. Mixed-metal target entries
  retain typed relations to out-of-domain centers.
- Clarified that 602,116 is a predefined three-dimensional d-block target
  corpus, not an all-metal CSD count.

## SM3 canonicalization replacement

- Added `attachment_set_key`, `donor_attachment_keys`, and
  `donor_rank_orbits` to the internal ligand object.
- Replaced source-index-dependent donor ordering with attachment-aware donor
  signatures containing donor element, source-index-free attachment identity,
  quantized metal–donor distance, and the sorted quantized donor–donor distance
  row.
- Defined legal ligand-instance and donor-rank orbits, Cartesian enumeration of
  residual legal permutations, lexicographic minimum-whole-record selection,
  candidate-count gating, and fail-closed missing-metadata behavior.
- Replaced the old small tie-test narrative with the locked N = 1,000,
  K = 100 three-arm challenge (300,000/300,000 exact matches) and the separate
  legal-orbit challenge (64,834/100,000 signature-only versus
  100,000/100,000 full exact).
- Added the nonautomorphic attachment negative control and the idealized
  fac/mer-[Co(NH3)3Cl3] relation-specificity control.

## Identity wording

- Restricted nested L0–L3 equality to validated, connectivity-supported
  records for which all four keys are formally defined.
- Replaced unconditional “exact snapshot identity” language with “finest
  encoded coordination-state equality at the stated serialization precision.”
- Described periodic MID/SID as a separate comparison axis, not another level
  of the molecular L0–L3 hierarchy.

## Reproducibility and table replacements

- Corrected SM8.5 references from nonexistent Tables S16–S18 and Figure S3 to
  the active Tables S7, S8, S2A/S2B, S10, and S14A/S14B.
- Replaced Supplementary Table S2 with S2A (general invariance, legal-orbit,
  attachment and fac/mer specificity controls) and S2B (synthetic/expanded
  record regressions).
- Split the census presentation into Supplementary Tables S9A and S9B; added
  oversized collective-π, confirmed-site overlap, unresolved periodic-edge,
  and edge-less audit-only rows.
- Corrected Table S11 so L0 is not described as a lossless Cartesian identity.
- Rewrote Tables S12–S13 to describe local orientation parity without assigning
  macroscopic chirality, enantiopurity, or a continuous chirality magnitude.
- Corrected the MID/SID sum to `12,965 + 2 × 23,587 = 60,139` and corrected the
  determinant-negative benchmark so it validates MID only.

## New Supplementary Figure S2

**Supplementary Figure S2. Search complexity and computational cost of
attachment-aware exact canonicalization.** Panel A gives the exact
legal-candidate-count distribution among 5,173 challenge-eligible graphs.
Panel B reports descriptive per-call median, 95th-percentile, and maximum wall
times for signature-only sorting and the full exact search. The timing is
explicitly presented as an implementation diagnostic, not a cross-machine
benchmark.

The complete S2 package is in `Supplementary_Figure_S2/` and includes SVG,
PDF, PNG, TIFF, caption, source CSV, generator, QA, and checksums.

## Rendering verification

The updated DOCX was exported through Microsoft Word to a 65-page PDF and all
65 pages were rendered for visual inspection. The added S2 artwork occupies a
dedicated page; table rows are prevented from splitting across pages and table
headers repeat after page breaks.
