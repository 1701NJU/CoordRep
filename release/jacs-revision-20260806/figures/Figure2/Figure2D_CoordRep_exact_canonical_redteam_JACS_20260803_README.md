# Figure 2d replacement: exact canonical labeling

## Purpose

This panel replaces the old hidden-duplicate count in Figure 2d. It shows the
specific symmetry red-team requested for the revised canonical claim: a single
homogeneous Cu4 cycle is supplied under all 24 metal input orders. The former
WL-only beta implementation gives three strings, whereas exact canonical
labeling collapses 24/24 orders to one CoordRep-State.

The bottom strip reports the broader nuisance-perturbation audits:

- 55 multinuclear record graphs x 100 variants = 5,500/5,500 passed.
- 33 haptic record graphs x 100 variants = 3,300/3,300 passed.

## Files

- `Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803.svg`: editable
  native vector; text remains text.
- `Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803.pdf`: publication
  vector with embedded TrueType fonts.
- `Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803.png`: 600 dpi
  preview.
- `Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803_source.csv`:
  machine-readable plotted values.
- `Figure2D_CoordRep_exact_canonical_redteam_JACS_20260803_source.py`:
  deterministic Matplotlib source.

## Data provenance

Values are copied from:

`revision_experiments/coordrep_v2beta_full_canonical_upgrade_20260803/revision_results/canonical_invariance_upgrade/canonical_invariance_audit_summary.json`

## Recommended caption sentence

**(d)** Exact canonical labeling removes input-order ambiguity in symmetric
expanded records. For a homogeneous Cu4 cycle, the preliminary WL-only beta
implementation produced three strings across the 24 possible metal orders,
whereas the exact procedure maps all 24/24 orders to one CoordRep-State. The
same invariance tests passed for 5,500/5,500 nuisance variants of 55 curated
multinuclear records and 3,300/3,300 variants of 33 curated haptic records.
