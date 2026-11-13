# Figure 4 artwork and source statistics

The complete Figure 4 preserves panels A/B of the submitted artwork. Only
C/D were redrawn with the checkpoint-backed November property statistics.
The point estimates, asymmetric 95% bootstrap endpoints and direction
concordances are in `Figure4CD_source_data.csv`.

- C: 34 graph-identical cis/trans pairs, dipole-difference MAE.
- D: 44 graph-identical different-shape pairs in 31 graph clusters, gap-difference
  MAE; the subset is conditioned on an observed gap difference above 0.10 eV.
- The confidence limits use 20,000 attributed-graph-cluster bootstrap resamples
  conditional on each frozen three-seed ensemble.
- Both plotted point values and all 20 MAE confidence endpoints were checked
  independently against the saved-weight statistics. Direction intervals are
  included in the CSV, not plotted as MAE intervals.

Use `Figure4_updated_20261114.png` or `.tif` for the full 600-dpi figure.
`Figure4_updated_native.png` retains the original 2704 by 2034 pixel canvas;
its A/B rectangle is pixel-identical to `Figure4_AB_original.png`.
The complete SVG/PDF preserves A/B as an embedded original raster; C/D remain
vector artwork. The independent C/D SVG and PDF files contain no raster images.
No native CSD coordinates or editable PyMOL scenes are included.

## Regenerate the artwork

Install Matplotlib, NumPy, Pillow and pypdf. Arial must be installed or supplied
from a licensed local font directory; font binaries are not redistributed.

```shell
python render_figure4.py --output rendered
python export_pdfs.py --input rendered
```

Optional: `--font-dir /path/to/Arial`, `--source source.csv`, or
`--ab-image Figure4_AB_original.png`. The default input files are next to the
script, so regeneration needs neither the manuscript DOCX nor CSD access.
The script redraws publication statistics; it does not run property models or
recompute bootstrap intervals. Model inference and scoring are documented in
the [November property release](https://github.com/1701NJU/CoordRep/releases/tag/coordrep-property-verification-20261111).

The August C/D artwork remains in the August package as historical evidence;
use the files here for the checkpoint-backed manuscript numbers. The original
A/B chemical targets and 31-pair chemical statistics are unchanged.
