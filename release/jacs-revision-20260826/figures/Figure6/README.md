# Figure 6 — molecular resolution and periodic local orientation

The current master is
`Figure6_Multiresolution_Periodic_Stereo_20260823.svg`. The PDF, 600 dpi
PNG/TIF, and preview were regenerated on 2026-08-28 with the established Figure
6 plotting logic after replacing the molecular panel-B aggregates with the
CoordRep 1.1.2rc3 rerun values. Panels D–E retain the frozen periodic-v3 data
and construction.

## Molecular rc3 boundary (panels A–C)

The frozen rc2 selection comprised 101,878 strict molecular records and 84,453
connectivity-supported records. The end-to-end April 2025 CSD rerun under rc3
yielded 101,794 and 84,449 records, respectively. These are eligibility
denominators and are not substituted silently for the frozen selection.

The active family cohort is exactly unchanged: 3,491 multi-record base-refcode
families contain the same 9,056 records under rc2 selection and rc3
re-encoding, with no missing L0–L3 identity key. Panel B therefore reports the
rc3 re-encoded values for that fixed cohort. Mean pairwise within-family
agreement is 11.2692%, 83.4860%, 96.5993%, and 97.3422% at L0–L3. Complete
concordance occurs for 357, 2,836, 3,350, and 3,378 families, respectively.
The five-member XEYVEC example retains the same 5 → 2 → 2 → 1 identity pattern.

The permission-safe aggregate comparison, method, and QC lock are under
`validation/molecular_rc3/`. Licensed row-level CSD mappings are not included.

## Periodic-v3 boundary (panels D–E)

Panels D–E use the separate frozen periodic-v3 MID/SID construction: 36,552
MIDs resolve to 60,139 SIDs, with 23,587 MIDs represented by exactly two local-
orientation SIDs. Within-entry pairing occurs in 7,721/9,189 entries not
flagged as Sohncke and 0/1,759 entries flagged as Sohncke. These values and the
400-entry/5,920-site re-expression audit were not changed by the molecular rc3
update.

`Figure6_manuscript_caption_20260828.txt` is the publication-facing caption;
`Figure6_caption_20260823.txt` retains the complete definitions, rc2-selection
versus rc3-eligibility boundary, and numerical detail. Licensed coordinate
extracts used for the displayed first-sphere renderings are not redistributed.
