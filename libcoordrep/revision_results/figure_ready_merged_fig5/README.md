# Figure 5 — Figure-Ready Data (Merged)

**Figure 5.** Full-CSD CoordRep-PathFinder maps coordination-geometry continua
and record-level structural diversity.

## Panels

| Panel | Title | Key file |
|-------|-------|----------|
| A | CoordRep-PathFinder workflow | `fig5A_workflow_summary.json` |
| B | Full-CSD scope audit | `fig5B_full_csd_scope_audit.csv` |
| C | CN5 TBPY–SPY geometry ridge | `fig5C_cn5_tbpy_spy_plot.csv` |
| D | Family-level structural polymorphism | `fig5D_family_polymorphism_summary.csv` |
| E | Casebook mini-cards | `fig5E_casebook_minicards.csv` |
| F | Quantitative output summary table | `fig5F_quantitative_summary_table.csv` |

## Data provenance

1. **All Fig. 5 numbers come from the full-CSD run** (1,413,222 entries scanned).
   The old 200K sample is no longer used.
2. **Fig. 5B** rejection categories come from the full-CSD filtering waterfall
   (`full_csd_scan_summary.json → waterfall`).
3. **Fig. 5C** CN5 atlas records: **14,897**.
4. **Fig. 5D** family counts:
   - 6,358 nontrivial L3 families
   - 3,739 multiple L1 families
   - 1,810 crossing boundary families
5. **Boundary records**: 17,791 (14.3%).
6. **No CSD raw coordinates / CIF / XYZ exported.**
7. Output files contain only refcodes, hashed keys, aggregate statistics,
   CShM values, and field summaries.

## Key numbers (for caption cross-check)

| Metric | Value |
|--------|-------|
| CSD release | 2024.3 |
| Total entries scanned | 1,413,222 |
| Valid CoordRep v1 records | 124,837 |
| Retention rate | 8.83% |
| CN4 atlas | 47,187 |
| CN5 atlas | 14,897 |
| CN6 atlas | 41,381 |
| CN5 intermediate (delta<1) | 42.2% |
| CN5 enrichment vs random | 4.4× |
| Nontrivial L3 families | 6,358 |
| Multi-L1 families | 3,739 |
| Boundary-crossing families | 1,810 |
| Boundary records | 17,791 (14.3%) |
| Top trajectory | Pd/CN4 LINMOL, span = 10.10 |

## License

CSD-derived analysis results only; raw coordinates are not redistributed.

## Files

```
figure_ready_merged_fig5/
├── README.md
├── fig5A_workflow_summary.json
├── fig5B_full_csd_scope_audit.csv
├── fig5B_full_csd_scope_audit.json
├── fig5B_outside_scope_categories.csv
├── fig5C_cn5_tbpy_spy_plot.csv
├── fig5C_cn5_ridge_annotations.json
├── fig5C_cn5_ridge_enrichment.csv
├── fig5D_family_polymorphism_summary.csv
├── fig5D_top_trajectory_case.json
├── fig5E_casebook_minicards.csv
├── fig5E_casebook_minicards.json
├── fig5F_quantitative_summary_table.csv
├── fig5_all_caption_numbers.json
├── fig5B_scope_audit_preview.png          (optional)
├── fig5C_cn5_tbpy_spy_hexbin_preview.png  (optional)
├── fig5D_family_polymorphism_bar_preview.png (optional)
└── fig5F_summary_table_preview.png        (optional)
```
