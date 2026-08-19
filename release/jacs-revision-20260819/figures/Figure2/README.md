# Corrected Figure 2 evidence package

This package replaces the old Figure 2D CShM atlas and the values formerly reported for EBAGAR/EBAGEV.

## What changed

- The old 4,449-record Ru/CN=6 atlas was retired because it used the superseded CShM objective.
- The background was rebuilt from the four completed shards in `corrected_csd_strict_full_scan_20260723_v2_hapticityfix`.
- The corrected strict cohort contains **4,365** Ru-centered CN=6 records; **4,334 (99.3%)** fall in the plotted window.
- EBAGAR: `S(Oh)=1.0148562824`, `S(TPr)=14.3127073886`.
- EBAGEV: `S(Oh)=0.9115136540`, `S(TPr)=13.5819441191`.
- The fac/mer claim is supported by the explicit trans-donor relation inventory, not by a large separation in CShM space.

## Scope and provenance

The atlas is the strict mononuclear, nondisordered, nonpolymeric, atom-resolved eta1 core emitted by CoordRep 1.1.2rc2 from the April 2025 CSD scan. CShM values use the corrected exact coordinate-space calculation. The row-level CSD-derived files are licensed internal derivatives and are not copied into this package. `*_Ru_CN6_hexbin.csv` contains only aggregate bin counts.

## Rebuild

Run `python generate_figure2.py` from this directory after supplying the licensed
strict-scan inputs referenced by the script. The script checks the expected
cohort size and both highlighted CShM vectors before writing source renders.

## Outputs

- `Figure2_Canonicalization_FINAL.svg`: manuscript vector artwork.
- `Figure2_Canonicalization_source_render.pdf`: reproducible source render.
- `Figure2_Canonicalization_preview_600dpi.png`: layout preview.
- `Figure2_source_data.csv`: compact panel data and exact highlighted values.
- `Figure2_Ru_CN6_hexbin.csv`: aggregate corrected atlas bins; no row-level CSD records.
- `Figure2_caption.txt`: current caption.
