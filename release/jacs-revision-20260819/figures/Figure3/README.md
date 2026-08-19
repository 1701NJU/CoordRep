# Corrected Figure 3: continuous shape-coordinate atlas

This folder contains the manuscript-ready Figure 3 rebuilt from the corrected coordinate-space CShM implementation. It deliberately replaces every number inherited from the earlier 65,017-record draft.

## Locked corpus and scope

- Fixed input cohort: 48,057 tmQM/GFN2-xTB records.
- Supported CN2-CN6 shape records: 35,485.
- Records plotted in panels A-C (CN4-CN6): 33,863.
- CN4: 16,682 total; 13,011 nearest SP and 3,671 nearest Td.
- CN5: 6,356 total; 5,319 nearest SPY and 1,037 nearest TBP.
- CN6: 10,825 total; 5,580 nearest Oh and 5,245 nearest TPr.
- Boundary rule: the unrounded top-two margin `Delta S = S(second) - S(best) < 1`.
- Boundary counts: CN4 39/16,682 (0.23%); CN5 693/6,356 (10.90%); CN6 294/10,825 (2.72%); combined 1,026/33,863 (3.03%).

## Provenance

- Primary source: `revision_experiments/results/corrected_cshm_xtb_atlas_v1/corrected_shape_atlas.csv`.
- Primary-source SHA-256: `76f1a7e77d73d0ae614b0b1781d1aa54236c4a4d622e286af59abfeb7d429348`.
- Source summary: `revision_experiments/results/corrected_cshm_xtb_atlas_v1/SUMMARY.json`.
- Source-summary SHA-256: `0040f987c1c04c6024a7e095efc44d9b0be6127cdc96d2a8a02ee0a1da3b5400`.
- Fixed-cohort input SHA-256 recorded by the source summary: `47b1882ecdb480cb974259897e0fa80dd67fbaef73c649e4a1c1a9d122f41a5a`.
- Figure generator: `generate_figure3.py` in this directory.

`Figure3_panel_source.csv` is the plot-ready long table.
`Figure3_summary_and_QC.csv` contains the plotted population and boundary
counts plus an explicitly internal frozen-to-corrected label-retention check.
The retention fields are not plotted and are not a chemical claim; they
document implementation sensitivity during correction. Release-level hashes
are recorded in `../../RELEASE_SHA256SUMS.txt`.

## Interpretation boundary

A nearest-reference label reports only which implemented ideal shape has the smaller corrected CShM for a fixed coordination number. The near-boundary descriptor is numerical and descriptive. Neither result is evidence for a reaction coordinate, a thermally driven transformation, or fluxional dynamics. No property table was loaded for this atlas.
